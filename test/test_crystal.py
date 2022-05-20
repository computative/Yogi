from numpy import array as ar
import sys, os, unittest
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../src")
from crystal import Crystal


class crystal_test(unittest.TestCase):

    # estend with other structures. 
    # Hard-code the desired structures.
    def setUp(self):
        self.type = "diamond"
        self.spp= ["Si","Ge"]
        self.atoms = { 
            "Si": [
                ar([0.00,0.00,0.00]),
                ar([0.00,5.5373/2,5.5373/2]),
                ar([5.5373/2,0.00,5.5373/2]),
                ar([5.5373/2,5.5373/2,0.00])
            ], "Ge": [
                ar([5.5373/4,5.5373/4,5.5373/4]),
                ar([5.5373/4,3*5.5373/4,3*5.5373/4]),
                ar([3*5.5373/4,5.5373/4,3*5.5373/4]),
                ar([3*5.5373/4,3*5.5373/4,5.5373/4])
            ]
        }

    # tests from_struct, from_coords, and the structures all in one go
    def test_from_methods(self):

        crl = Crystal().from_struct(
            { "type": self.type, "spp": self.spp }
        )

        crl_atoms = {}
        for sp_test, sp_std in zip(crl.nuclei, self.atoms):
            items_std = [tuple(item) for item in self.atoms[sp_std]]
            items_test = [tuple(item) for item in crl.nuclei[sp_test].values()]
            self.assertEqual( sp_test, sp_std )
            self.assertEqual( set( items_std ), set( items_test ) )



if __name__ == '__main__':
    unittest.main()