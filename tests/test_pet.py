from sbtn_leaf.PET import calculate_PET_location_based


def test_january_daylight_duration_does_not_raise_index_error():
    """Ensure the daylight lookup uses a valid month index for January."""

    temps = [10.0] * 12

    # A successful calculation implies no IndexError was raised for January.
    pet_values = calculate_PET_location_based(temps, 2024, 0.0)

    assert len(pet_values) == 12
