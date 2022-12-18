from numpy import array as ar

class Structure:
    a = {
        "Si" : {"Si": 5.431, "Ge": 5.5373},
        "Ge" : {"Si": 5.431}
    }

    atoms = {
        "sp": [[]]
    }
    rep = [
        [1,0,0],
        [0,1,0],
        [0,0,1]
    ]
