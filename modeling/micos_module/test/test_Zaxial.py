import IMP
import IMP.core
import IMP.algebra
import IMP.micos
import IMP.test
import IMP.atom

def setup_p(m, mass,name,vec):
    p = m.add_particle(str(name))
    p = IMP.core.XYZ.setup_particle(m, p)
    p.set_coordinates(IMP.algebra.Vector3D(*vec))
    p = IMP.atom.Mass.setup_particle(m,p,mass)
    return p

class Tests(IMP.test.TestCase):
    """ Test ZAxial Restraint"""

    def setUp(self):
        IMP.test.TestCase.setUp(self)
        IMP.set_log_level(IMP.SILENT)
        IMP.set_check_level(IMP.NONE)
        self.m = IMP.Model()
        self.particle_coordinates = [
            (-10,5,-100),
            (0, 0, 10),
            (0,3,0)
        ]
        self.mass = [5.0, 10.0, 2.0]
        self.om = 0
        self.cj = 3
        self.sigma = 1
        self.answers = 7.04

        self.particles = [setup_p(self.m, self.mass[i],i, self.particle_coordinates[i]) for i in
                          range(len(self.particle_coordinates))]

    def test_zaxial_restraint(self):
        """Test restraint"""
        res = IMP.micos.ZAxialRestraint(self.particles, self.cj, self.om,self.sigma, 'average')
        val = res.unprotected_evaluate(None)
        self.assertAlmostEqual(val, self.answers, places = 0) # this will check your ans with the computed score

if __name__ == '__main__':
    IMP.test.main()
