from dataclasses import dataclass, field
from gas import Gas
# from .strip import Strip


@dataclass
class Chamber:
    name : str
    nickname : str
    strips : dict = field(init=False, default_factory=dict)

    def create_layers(self, config):
        self.strips = {}
        for layer in str(config[self.name]['layers']).replace(" ", "").split(','):
            self.strips[layer] = {} 
            for readout in str(config[self.name]['readouts']).replace(" ", "").split(','):
                self.strips[layer][readout] = []


def create_chambers(config):
    keys = list(config.keys())
    i = 0
    chambers = []
    while keys[i] != "general":
        chambers.append(Chamber(keys[i], config[keys[i]]["name"]))
        chambers[i].create_layers(config)
        print(f"- {chambers[i]}")
        i += 1
    return chambers
    
    
def main() -> None:
    chamber = Chamber("ExMe")
    print(chamber)


if __name__ == "__main__":
    main()

