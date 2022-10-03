"""
A hit class
"""


class Hit:
    def __init__(self,coord,is_sidex,timestamp,timestamp_event,value):
        self.coord = coord
        self.is_sidex = is_sidex
        self.timestamp = timestamp
        self.timestamp_event = timestamp_event
        self.value = value
