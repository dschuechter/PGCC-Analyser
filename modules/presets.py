def import_preset(name):
    preset = [0]
    if name == "testing":
        preset = ["testing",0,30,-30,30,False,0,10,False, True, True, 0.7]
    if name == "fullsky":
        preset = ["fullsky",0,360,-90,90,False,0,10,False, True, True, 0.7]
    if name == "orion_nebula":
        preset = ["orion_nebula",200,215,-27,-15,False,0,0, False,True, True, 0.7]
    if name == "draco":
        preset = ["draco",84.1,95.9,32.9,42.8,False,0,0,False,True,True,0.7]
    if name == "magellanic_clouds":
        preset = ["magellanic_clouds",270,310,-50,-25,False,0,0,False, True,True, 0.7]
    if name == "SMC":
        preset = ["SMC",300.29,305.05,-46.86,-42.08,False,0,0,False, True,True,0.7]
    if name == "LMC":
        preset = ["LMC",273.94,283.12,-37.69,-29.51,False,0,0,False, True,True,0.7]
    return preset
