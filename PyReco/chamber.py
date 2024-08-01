from dataclasses import dataclass
from gas import Gas
# from .strip import Strip


@dataclass
class Chamber:
    name : str
    # gas : Gas

def create_chambers(config):
    keys = list(config.keys())
    i = 0
    chambers = []
    while keys[i] != "general":
        chambers.append(Chamber(config[keys[i]]["name"]))
        print(f"- {chambers[i]}")
        i += 1
    return chambers
    
    

def main() -> None:
    # gas = Gas(name = "Ar:Cf4:Iso", fraction="88:10:2", drift_velocity=0.105)
    # chamber = Chamber("ExMe", gas)
    chamber = Chamber("ExMe")
    print(chamber)


if __name__ == "__main__":
    main()

