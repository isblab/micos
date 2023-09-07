import IMP
import IMP.core
import IMP.algebra
import IMP.micos
import IMP.test

def setup_p(m, name, vec):
    p = m.add_particle(str(name))
    p = IMP.core.XYZ.setup_particle(m, p)
    p.set_coordinates(IMP.algebra.Vector3D(*vec))
    return p

class Tests(IMP.test.TestCase):
    """ Test ZAxial Restraint"""

    def setUp(self):
        IMP.test.TestCase.setUp(self)
        IMP.set_log_level(IMP.SILENT)
        IMP.set_check_level(IMP.NONE)
        self.m = IMP.Model()
        self.particle_coordinates = [
            (-10,5,-10.09),
            (0, 0, 10.4),
            (0,3,0.33),
            (1,2,-50),
            (2,3,46)
        ]

        self.lower_bound = -.4
        self.upper_bound = 0
        self.sigma = 1
        self.answers = 32

        self.particles = [setup_p(self.m, i, self.particle_coordinates[i]) for i in
                          range(len(self.particle_coordinates))]

    def test_zaxial_restraint(self):
        """Test restraint"""
        res = IMP.micos.ZAxialProteinwiseRestraint(self.particles, self.lower_bound,self.upper_bound,self.sigma)
        val = res.unprotected_evaluate()
        self.assertAlmostEqual(val, self.answers, places = 0) # this will check your ans with the computed score

if __name__ == '__main__':
    IMP.test.main()
