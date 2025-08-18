import IMP
import IMP.core
import IMP.algebra
import IMP.micos
import IMP.test

def setup_p(m, name, vec):
    p = m.add_particle(str(name))
    p = IMP.core.RigidBody.setup_particle(m, p, IMP.algebra.ReferenceFrame3D())
    p.set_coordinates(IMP.algebra.Vector3D(*vec))
    return p

class Tests(IMP.test.TestCase):
    """ Test surface localization Restraint"""

    def setUp(self):
        IMP.test.TestCase.setUp(self)
        IMP.set_log_level(IMP.SILENT)
        IMP.set_check_level(IMP.NONE)
        self.m = IMP.Model()
        self.particle_coordinates = [
            (-10,5,0),
            (0, 0, 1),
            (0,3,0)
        ]

        self.r = 7
        self.sigma = 1
        self.max_limit = 5
        self.answers = 21.47

        self.particles = [setup_p(self.m, i, self.particle_coordinates[i]) for i in
                          range(len(self.particle_coordinates))]

    def test_surface_localization(self):
        """Test restraint"""
        res = IMP.micos.SurfaceLocalizationRestraint(self.particles, self.r, self.sigma, self.max_limit)
        val = res.unprotected_evaluate(None)
        self.assertAlmostEqual(val, self.answers, places = 0) # this will check your ans with the computed score

if __name__ == '__main__':
    IMP.test.main()
