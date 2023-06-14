import IMP
import IMP.core
import IMP.algebra
import IMP.micos
import IMP.test

def setup_p(m, name, vec):
    p = m.add_particle(str(name))
    p = IMP.core.XYZ.setup_particle(m, p)
    p.set_coordinates(IMP.algebra.Vector3D(*vec))
    print( IMP.core.XYZ(p).get_coordinate(2))
    return p

class Tests(IMP.test.TestCase):
    """ Test Z axial Restraint"""

    def setUp(self):

        IMP.test.TestCase.setUp(self)
        IMP.set_log_level(IMP.SILENT)
        IMP.set_check_level(IMP.NONE)
        self.m = IMP.Model()
        self.particle_coordinates = [
            (-10,5,0),
            (0, 0, 15),
            (0,3,-4),
            (4,2,-60)
        ]

        self.sigma = 1
        self.answers = 65.0
        self.particles = [setup_p(self.m, i, self.particle_coordinates[i]) for i in
                          range(len(self.particle_coordinates))]

    def zaxial(self):
        """Test z axial restraint"""

        res = IMP.micos.ZAxialRestraint(self.particles,self.sigma)
        val = res.unprotected_evaluate(None)
        self.assertAlmostEqual(val, self.answers, places = 0) # this will check your ans with the computed score

if __name__ == '__main__':
    IMP.test.main()
