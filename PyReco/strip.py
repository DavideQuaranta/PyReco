from dataclasses import dataclass


@dataclass
class Strip:
    """Represents a strip with its ID, charge, and time.

    Attributes:
        id (int): The ID of the strip.
        charge (float): The charge of the strip.
        time (float): The time associated with the strip.
    """

    id: int
    charge: float
    charge_error: float
    time: float
    time_error: float


def main() -> None:
    strip = Strip(id = 100, charge = 100, charge_error = 20, time = 250, time_error = 25)
    print(strip)


if __name__ == "__main__":
    main()
