import numpy as np

class Cluster:
    """Represents a cluster of strips associated with a chamber.

    Attributes:
        chamber (Chamber): The chamber associated with the cluster.
        strips (numpy.ndarray): The array of strips in the cluster.
        size (int): The size of the cluster.
        total_charge (float): The total charge of the cluster.
    """

    def __init__(self, chamber):
        """Initializes a Cluster object with the given chamber.

        Args:
            chamber (Chamber): The chamber associated with the cluster.
        """
        self.chamber = chamber
        self.strips = np.empty(0)  # Initialize an empty NumPy array
        self.size = 0
        self.total_charge = 0
    
    def AddStrip(self, strip):
        """Adds a strip to the cluster.

        Args:
            strip (Strip): The strip to add to the cluster.
        """
        self.strips = np.append(self.strips, strip)  # Append strip to the NumPy array
        self.size += 1
        self.total_charge += strip.charge

    def __add__(self, other):
        """Concatenates two clusters if they belong to the same chamber.

        Args:
            other (Cluster): The cluster to concatenate with.

        Returns:
            Cluster: The concatenated cluster if successful, NotImplemented otherwise.
        """
        if isinstance(other, Cluster) and self.chamber == other.chamber:
            new_cluster = Cluster(self.chamber)
            new_cluster.strips = np.concatenate((self.strips, other.strips))  # Concatenate NumPy arrays
            new_cluster.size = self.size + other.size
            new_cluster.total_charge = self.total_charge + other.total_charge
            return new_cluster
        return NotImplemented
    
    def isGoodCluster(self):
        """Checks if the cluster is considered good based on chamber parameters.

        Returns:
            bool: True if the cluster is good, False otherwise.
        """
        min_strips = self.chamber.min_strips_in_cluster
        max_strips = self.chamber.max_strips_in_cluster
        min_charge = self.chamber.min_cluster_charge
        
        return min_strips <= self.size <= max_strips and self.total_charge >= min_charge
        
    def Center(self):
        """Returns the center of the cluster."""
        return np.mean([strip.id for strip in self.strips])
    
    def GetDeltaT(self):
        """Returns the time difference between the earliest and latest time stamps of the strips."""
        time_stamps = np.array([strip.time for strip in self.strips])
        return np.max(time_stamps) - np.min(time_stamps) + 0.00000000000001

    def Centroid(self):
        """Returns the centroid of the cluster."""
        pitch = self.chamber.pitch
        qtot = self.total_charge
        strip_ids = np.array([strip.id for strip in self.strips])
        strip_charges = np.array([strip.charge for strip in self.strips])
        centroid = np.sum(strip_charges * strip_ids * pitch)
        return centroid / qtot if qtot != 0 else 0  # Avoid division by zero   
    
    def Print(self):
        """Prints information about the cluster."""
        print(f"Cluster size = {self.size}")
        for strip in self.strips:
            strip.Print()
