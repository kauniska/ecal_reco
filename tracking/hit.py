"""
Class implementing the hits.
Contains the coordinates of the hit on the concerned side, a bool is_sidex True if the concerned side is on the x plane,
the timestamp of the hit, the global timestamp of the event and the value of the hit

"""
from track_reconstruction import is_sidex
from track_reconstruction import mapping_2D

class Hit:
    
    def __init__(self):
        """
        The construction is void, the parameters are set through set_parameters or set_data_frame
        """
        self.coord = None
        self.is_sidex = None
        self.timestamp = None
        self.timestamp_event = None
        self.value = None

    def set_parameters(self,coord,is_sidex,timestamp,timestamp_event,value):
        """
        Sets the parameters manually, the coordinates have to be already transformed
        in the final form though the function track_reconstruction.mapping_2D
        """
        self.coord = coord
        self.is_sidex = is_sidex
        self.timestamp = timestamp
        self.timestamp_event = timestamp_event
        self.value = value

    def set_data_frame(self,data_frame, local_index):
        """
        Sets the parameters with an event of the data frame and the index of the hit local to the given event
        """
        self.coord = mapping_2D(data_frame['tofpet_id'][local_index],data_frame['tofpet_channel'][local_index])
        if is_sidex(data_frame['tofpet_id'][local_index]):
            self.is_sidex = True
        else:
            self.is_sidex = False
        self.timestamp = data_frame['timestamp'][local_index]
        self.timestamp_event = data_frame['timestamp_event'][local_index]
        self.value = data_frame['value'][local_index]
        
        
        
