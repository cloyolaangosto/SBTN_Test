import datetime as dt

import polars as pl
import pytest

from sbtn_leaf.PET import (
    calculate_PET_crop_based,
    calculate_PET_location_based,
    create_KC_Curve,
)


def test_january_daylight_duration_does_not_raise_index_error():
    """Ensure the daylight lookup uses a valid month index for January."""

    temps = [10.0] * 12

    # A successful calculation implies no IndexError was raised for January.
    pet_values = calculate_PET_location_based(temps, 2024, 0.0)

    assert len(pet_values) == 12


def test_create_kc_curve_respects_stage_lengths_without_rollover():
    """The generated Kc curve should honour configured stage durations."""

    def _absolute_day_table():
        base = dt.date(2021, 1, 1)
        rows = []
        for offset in range(365):
            current = base + dt.timedelta(days=offset)
            rows.append(
                {
                    "Date": f"{current.day}-{current.strftime('%b')}",
                    "Day_Num": offset + 1,
                    "Day": current.day,
                    "Month": current.month,
                }
            )
        return pl.DataFrame(rows)

    abs_table = _absolute_day_table()
    planting_date = "15-Mar"
    planting_day = abs_table.filter(pl.col("Date") == planting_date).select("Day_Num").item()

    crop_table = pl.DataFrame(
        {
            "Climate_Zone": ["TestZone"],
            "Crop": ["TestCrop"],
            "K_ini": [0.5],
            "K_mid": [1.1],
            "K_Late": [0.8],
            "Initial_days": [2],
            "Dev_Days": [3],
            "Mid_Days": [4],
            "Late_days": [2],
            "Planting_Greenup_Date": [planting_date],
            "Soil_Cover_Period": [0],
            "SCP_Starts": [0],
            "SCP_End": [0],
        }
    )

    kc_curve = create_KC_Curve(
        "TestCrop",
        "TestZone",
        crop_table=crop_table,
        abs_date_table=abs_table,
    )

    stage_lengths = {
        1: 1,
        2: 2,
        3: 3,
        4: 4,
        5: 2,
    }

    for stage_id, expected in stage_lengths.items():
        assert (
            kc_curve.filter(pl.col("Stage_id") == stage_id).height == expected
        ), f"Stage {stage_id} duration mismatch"

    planting_stage_days = kc_curve.filter(pl.col("Stage_id") == 1)["day_num"].to_list()
    assert planting_stage_days == [planting_day]

    initial_stage_days = kc_curve.filter(pl.col("Stage_id") == 2)["day_num"].to_list()
    assert initial_stage_days[0] == planting_day + 1
    assert planting_day not in initial_stage_days

    active_days = kc_curve.filter(pl.col("Stage_id") > 0)["day_num"].to_list()
    assert max(active_days) == planting_day + sum(stage_lengths.values()) - 1
    assert max(active_days) <= 365


def test_calculate_pet_crop_based_respects_leap_year_days():
    """Daily PET totals for February should honour leap-year month lengths."""

    def _absolute_day_table():
        base = dt.date(2021, 1, 1)
        rows = []
        for offset in range(365):
            current = base + dt.timedelta(days=offset)
            rows.append(
                {
                    "Date": f"{current.day}-{current.strftime('%b')}",
                    "Day_Num": offset + 1,
                    "Day": current.day,
                    "Month": current.month,
                }
            )
        return pl.DataFrame(rows)

    crop_table = pl.DataFrame(
        {
            "Climate_Zone": ["TestZone"],
            "Crop": ["LeapCrop"],
            "K_ini": [29 / 28],
            "K_mid": [29 / 28],
            "K_Late": [29 / 28],
            "Initial_days": [364],
            "Dev_Days": [0],
            "Mid_Days": [0],
            "Late_days": [0],
            "Planting_Greenup_Date": ["1-Jan"],
            "Soil_Cover_Period": [0],
            "SCP_Starts": [0],
            "SCP_End": [0],
        }
    )

    abs_table = _absolute_day_table()
    monthly_temps = [10.0] * 12

    results = calculate_PET_crop_based(
        "LeapCrop",
        "TestZone",
        monthly_temps,
        2024,
        0.0,
        crop_table=crop_table,
        abs_date_table=abs_table,
    )

    assert "PET_Annual" in results

    feb_total = (
        results["PET_Daily"]
        .filter(pl.col("Month") == 2)
        .select(pl.col("PET_Daily").sum())
        .item()
    )

    feb_monthly_input = calculate_PET_location_based(monthly_temps, 2024, 0.0)[1]

    assert feb_total == pytest.approx(feb_monthly_input, rel=1e-9)

    monthly_total = results["PET_Monthly"].select(pl.col("PET_Monthly").sum()).item()

    assert results["PET_Annual"] == pytest.approx(monthly_total, rel=1e-9)
