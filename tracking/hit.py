"""
Class implementing the hits.
Contains the coordinates of the hit on the concerned side, a bool is_sidex True if the concerned side is on the x plane,
the timestamp of the hit, the global timestamp of the event and the value of the hit

"""
from track_reconstruction import is_sidex
from track_reconstruction import mapping_2D

class Hit:

    def __init__(self,*args):
        """Creates a Hit from arguments.
            If zero argument, then it's an empty instance
            If two arguments (event from data frame and local index), then it extracts the parameters from the data frame
            If five arguments (coordinates, is_sidex, timestamp, timestamp_event, value), then it simply writes all the parameters
            
        Args:
            coordinates : 2D coordinates of the hit on the concerned side
            is_sidex : bool = True if the concerned side is on the x plane
            timestamp : timestamp of the hit
            timestamp_event : global timestamp of the hole event
            value : energy of the hit
        """
        if len(args) == 0:
            self.coord = None
            self.is_sidex = None
            self.timestamp = None
            self.timestamp_event = None
            self.value = None
        elif len(args) == 2:
            self.coord = mapping_2D(args[0]['tofpet_id'][args[1]],args[0]['tofpet_channel'][args[1]])
            if is_sidex(args[0]['tofpet_id'][args[1]]):
                self.is_sidex = True
            else:
                self.is_sidex = False
            self.timestamp = args[0]['timestamp'][args[1]]
            self.timestamp_event = args[0]['timestamp_event']
            self.value = args[0]['value'][args[1]]
        elif len(args) == 5:
            self.coord = args[0]
            self.is_sidex = args[1]
            self.timestamp = args[2]
            self.timestamp_event = args[3]
            self.value = args[4]
        else:
            raise ValueError("not the correct number of arguments given")

    def print(self):
        if self.is_sidex:
            plane = "x"
        else:
            plane = "y"
        print("Coordinates : " , self.coord, "on", plane , "plane")
        print("Timestamp : ", self.timestamp)
        print("Timestamp event : ", self.timestamp_event)
        print("Value : ", self.value)
        
        

