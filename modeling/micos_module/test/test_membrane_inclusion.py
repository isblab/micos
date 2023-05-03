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
    """ Test Membrane Inclusion Restraint"""

    def setUp(self):
        IMP.test.TestCase.setUp(self)
        IMP.set_log_level(IMP.SILENT)
        IMP.set_check_level(IMP.NONE)
        self.m = IMP.Model()
        self.particle_coordinates = [
            (-10, 5, 0),
            (0, 0, 1),
            (0, 3, 0),
            (5, 0, 0),
            (10, 8, 7),
            (3, 6, 9),
            (-4, -2.5, -4)
        ]

        self.r = 7
        self.R = 14
        self.sigma = 1

        self.particles = [setup_p(self.m, i, self.particle_coordinates[i]) for i in
                          range(len(self.particle_coordinates))]

    def test_membrane_inclusion(self):
        """Test restraint"""
        res = IMP.micos.MembraneInclusionRestraint(self.particles, self.sigma, self.R, self.r)
        val = res.unprotected_evaluate(None)


if __name__ == '__main__':
    IMP.test.main()
