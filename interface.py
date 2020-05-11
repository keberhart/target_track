#!/usr/bin/env python3

from consolemenu import SelectionMenu
from antenna import Antenna
from target import Target
import sys

class Interface():

    def __init__(self):
        self.system = Antenna()
        ant_list = self.system.get_list()
        self.target = Target()
        target_list = self.target.get_list()

        title = "INP Generator for Orbit ACUs"
        subtitle = "Choose your Antenna"
        # A menu to select our antenna
        selection = SelectionMenu.get_selection(ant_list, title=title, subtitle=subtitle)
        if selection == len(ant_list):
            sys.exit()

        sys_name = ant_list[selection]

        subtitle = "Using " + sys_name + ", Choose the target"
        # A menu to select our target
        selection = SelectionMenu.get_selection(target_list, title=title, subtitle=subtitle)
        if selection == len(target_list):
            sys.exit()

        target_name = target_list[selection]

        self.system.ready(sys_name)
        self.target.ready(target_name)
        self.target.ready_antenna(self.system.lat, self.system.lon, self.system.alt)
        self.target.generate_report()

        if self.target.rise_time is None:
            output = "\n\n  "+target_name+" doesn't rise in the next 12 hours...\n\n"
            print(output)
            sys.exit()
        else:
            output = "\n\n "+target_name+" will be up at "+self.target.rise_time+"\n"
            print(output)
            output = " The report ends or "+target_name+" fades at "+self.target.fade_time+"\n\n"
            print(output)
            filename = "/home/scheduler/Desktop/"+sys_name+"_"+target_name+".orb"
            output = " Saving INP file to "+filename+"\n\n"
            print(output)
            with open(filename, 'w') as f:
                for line in self.target.report:
                    f.write(line)
            sys.exit()


if __name__ == "__main__":
    runner = Interface()

