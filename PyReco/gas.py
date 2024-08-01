from dataclasses import dataclass, field


def initialize_dict(name, fraction):
    components = name.split(":")
    fractions = fraction.split(":")
    composition = {}

    for c,f in zip(components, fractions):
        composition[c] = f
    
    return composition


@dataclass(slots=True)
class Gas:
    """Represents the gas used in the Micromegas

    Attributes:
        name (str): The components of the mixture separated by colons.
        fraction (str): The mass fraction of the mixture components separated by colons.
        composition: dict[str, str]: A dictionary containing gas composition.
        drift_velocity (float): The drift velocity of the gas mixture in mm/ns.
    """
    name: str
    fraction: str
    drift_velocity: float
    composition: dict[str, str] = field(init = False, repr=False)

    def __post_init__(self):
        self.composition = initialize_dict(self.name, self.fraction)


def main() -> None:
    gas = Gas(name="Ar:CF4:Iso", fraction="88:10:2", drift_velocity=0.105)
    print(gas)
    print(gas.composition)


if __name__ == "__main__":
    main()