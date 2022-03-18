def DiamondCell(start_Id, species, end_plate):
        """DEFINE SPECIES AND COUNT OF ADDED NUCLEI"""
        (sp1, sp2) = species
        nuclei_added_to_cell = 8
        """ INSERT NUCLEI NOT IN END PLATE """
        coords = {
            sp1: {
                start_Id + 0: (),
                start_Id + 1: (),
                start_Id + 2: (),
                start_Id + 3: (),
            }, sp2: {
                start_Id + 4: (),
                start_Id + 5: (),
                start_Id + 6: (),
                start_Id + 7: ()
            }
        }
        """ INSERT END PLATES IF END BLOCK """
        for i in end_plate:
            if i == 0:
                coord1, coord2 = 1
            elif i == 1:
                coord1, coord2 = 1
            else:
                coord1, coord2 = 1
            coords.update({
                sp1: {
                    start_Id + 8: coord1,
                }, sp2: {
                    start_Id + 9: coord1
            }})
            nuclei_added_to_cell + 2
        return coords, nuclei_added_to_cell