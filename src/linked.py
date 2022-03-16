
# can i get rid of nuclei alltogether?
# for hamiltonian, I can loop over 
class linked:
    
    # sets a list with IDs and attributes
    # leave links empty
    def __init__(self, attribs):
        self.data = []
        for attr,_id in zip(attribs,range(len(attribs))):
            self.data.append([_id,[],attr])

    def __iter__(self):
        return iter(self.data)
    
    def __next__(self,itr):
        return next(itr)

    # link ID with targets
    def link(self,_id,targets):
        if _id in targets:
            raise ValueError("Attribute linked to itself.")
        self.data[_id][1] = targets

    #unlik all attributes
    def unlink(self,_id):
        self.data[_id][1] = []

    # link _id to attributes with values
    # in targets in column number col
    def link_attr(self,_id,col,targets):
        _list = []
        for entry in self.data:
            attr = entry[2][col]
            if attr in targets:
                __id = entry[0]
                # do not link an attr to itself
                if __id == _id:
                    raise ValueError("Attribute link error.")
                _list.append(__id)
        self.data[_id][1] = _list

    # return all occurences that has attribute value in column col
    def get_attr(self, col, value):
        temp = []
        for entry in self.data:
            if np.all(entry[2][col] == value):
                temp.append(entry)
        return temp

    # entries linked to _id
    def retrieve(self, _id):
        target_ids = self.data[_id][1]
        _list = []
        for __id in target_ids:
            _list.append(self.data[__id])
        return _list

    # return all ids that has attribute value in column col
    def search(self, col, value):
        temp = []
        for entry in self.data:
            if np.all(entry[2][col] == value):
                temp.append(entry[0])
        return temp
