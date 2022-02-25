from utils.polygonsphere import CoordsConverter

class PointList(object):
    
    def __init__(self):
        self.coo = []
        self.index = {}
        self.cnt
        
    def add_points(self,point_list):
        conv = CoordsConverter()
        for pnt in point_list:
            self.coo.append(pnt)
            xyz = conv.get_xyz(pnt)
            for key in self.index.keys():
                if key == 'll':
                    self.index['xy'].add(self.cnt,
                                         (pnt[0][0],pnt[0][1],
                                          pnt[0][0],pnt[0][1]))
                    self.cnt += 1
                elif key == 'xy':
                    self.index['xy'].add(self.cnt,
                                         (xyz[0][0],xyz[0][1],
                                          xyz[0][0],xyz[0][1]))
                    self.cnt += 1