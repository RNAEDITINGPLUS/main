# -* - coding: UTF-8 -* -  
import ConfigParser

conf = ConfigParser.ConfigParser()
conf.read("config.conf")

conf.set("splicing", "mesSCT", raw_input("mesSCT: "))
conf.set("splicing", "mesVT", raw_input("mesVT: "))
#conf.set("splicing", "mesIVT", raw_input("mesIVT: "))
#conf.set("splicing", "bpThreshold", raw_input("bpThreshold: "))
#conf.set("splicing", "bpVarThreshold", raw_input("bpVarThreshold: "))

conf.write(open("config.conf","w"))
