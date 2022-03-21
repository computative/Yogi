def DiamondCell(origin, start_Id, spp, end_plate):
        """DEFINE SPECIES AND COUNT OF ADDED NUCLEI"""
        nuclei_added_to_cell = 8
        x,y,z = origin
        sp1, sp2 = spp

        coords = {}
        coords.update({sp1: {}}); coords.update({sp2: {}})
        
        """ INSERT NUCLEI NOT IN END PLATE """

        coords[sp1].update({
                start_Id + 0: (x+0,y+0,z+0),
                start_Id + 4: (x+0,y+1/2,z+1/2),
                start_Id + 1: (x+1/2,y+0,z+1/2),
                start_Id + 5: (x+1/2,y+1/2,z+0)}) 
        coords[sp2].update( {
                start_Id + 2: (x+1/4,y+1/4,z+1/4),
                start_Id + 6: (x+1/4,y+3/4,z+3/4),
                start_Id + 3: (x+3/4,y+1/4,z+3/4),
                start_Id + 7: (x+3/4,y+3/4,z+1/4) } )
        """ INSERT END PLATES IF END BLOCK """
        
        if len(end_plate) == 1:
            nuclei_added_to_cell += 2
            if 0 in end_plate:
                coords[sp1].update({
                    start_Id + 8:(x+1,y+0,z+0),
                    start_Id + 9:(x+1,y+0.5,z+0.5) } )
            elif 1 in end_plate:
                coords[sp1].update({
                    start_Id + 8: (x+0,y+1,z+0),
                    start_Id + 9: (x+0.5,y+1,z+0.5) } )
            elif 2 in end_plate:
                coords[sp1].update({
                    start_Id + 8: (x+0,y+0,z+1),
                    start_Id + 9: (x+0.5,y+0.5,z+1) } )
        elif len(end_plate) == 2:
            nuclei_added_to_cell += 5
            if 0 in end_plate and 1 in end_plate:
                coords[sp1].update({
                    start_Id + 8:  (x+1,y+0,z+0),
                    start_Id + 9:  (x+0,y+1,z+0),
                    start_Id + 10: (x+1,y+1,z+0),
                    start_Id + 11: (x+1,y+0.5,z+0.5),
                    start_Id + 12: (x+0.5,y+1,z+0.5)} )
            elif 1 in end_plate and 2 in end_plate:
                coords[sp1].update({
                    start_Id + 8:  (x+0,y+1,z+0),
                    start_Id + 9:  (x+0,y+0,z+1),
                    start_Id + 10: (x+0,y+1,z+1),
                    start_Id + 11: (x+0.5,y+1,z+0.5),
                    start_Id + 12: (x+0.5,y+0.5,z+1) } )
            elif 0 in end_plate and 2 in end_plate:
                coords[sp1].update({
                    start_Id + 8: (x+1,y+0,z+0),
                    start_Id + 9: (x+0,y+0,z+1),
                    start_Id + 10:(x+1,y+0,z+1),
                    start_Id + 11:(x+1,y+0.5,z+0.5),
                    start_Id + 12:(x+0.5,y+0.5,z+1) } )
        elif len(end_plate) == 3:
            nuclei_added_to_cell + 10
            coords[sp1].update({ 
                start_Id + 8:  (x+1,y+0,z+0),
                start_Id + 9:  (x+0,y+1,z+0),
                start_Id + 10: (x+0,y+0,z+1),
                start_Id + 11: (x+0,y+1,z+1),
                start_Id + 12: (x+1,y+0,z+1),
                start_Id + 13: (x+1,y+1,z+0),
                start_Id + 14: (x+1,y+1,z+1),
                start_Id + 15: (x+1,y+0.5,z+0.5),
                start_Id + 16: (x+0.5,y+1,z+0.5),
                start_Id + 17: (x+0.5,y+0.5,z+1) } )

        return coords, nuclei_added_to_cell

def _diamond(self, dims, spp):
    
    # the method must contain empty dict and
    # nx,ny,nz number of blocks and start_Id

    sp1, sp2 = spp
    nuclei = {}
    nuclei[sp1] = {}
    nuclei[sp2] = {}
    nx,ny,nz = dims
    start_Id = 0

    for i in range( nx ):
        for j in range( ny ):
            for k in range( nz ):
                # end_plate sees
                end_plate =  np.argwhere(
                    np.array([nx - i,ny - j,nz - k]) == 1 ).T[0]
                origin = (i,j,k)
                cell, added = Cell(origin, start_Id,
                    species, end_plate)
                start_Id += added
                for sp in np.unique(np.array(
                        self.settings["spp"])):
                    nuclei[sp].update(cell[sp])
    return nuclei
