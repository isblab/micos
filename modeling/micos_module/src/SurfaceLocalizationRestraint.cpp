#include <math.h>
#include <IMP/core/XYZ.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/micos/SurfaceLocalizationRestraint.h>
#include "IMP/Model.h"

IMPMICOS_BEGIN_NAMESPACE

SurfaceLocalizationRestraint::SurfaceLocalizationRestraint(IMP::ParticlesTemp plist, double r, double sigma, double max_limit) :

		Restraint(plist[0]->get_model(), "SurfaceLocalizationRestraint %1%"),
		plist_(plist),
		rsq_(r * r),
		sigma_(sigma),
		max_limit_ (max_limit) {}


double SurfaceLocalizationRestraint::unprotected_evaluate(IMP::Model *m, IMP::Particle* p, IMP::DerivativeAccumulator* accum) const {

	// check if rigid body
   	IMP_USAGE_CHECK(core::RigidBody::get_is_setup(m, p),
							"Particle is not a rigid body");


		double x = IMP::core::XYZ(p).get_coordinate(0);  // The global coordinates of the rb
		double y = IMP::core::XYZ(p).get_coordinate(1);
		// Calculate the coordinate-wise distance from the center

		double radial = (x * x) + (y * y);  // radial = (euclidean distance from center) ^ 2

		if (accum){};

		if (sqrt(rsq_) - sqrt(radial) > max_limit_) {  // Continue only if the distance greater than the max value
				double deviation = rsq_ + radial - 2 * sqrt(radial * rsq_);
				// deviation = squared difference of distance of particle and the inner radius
				return fabs(deviation)/sigma_;
		}
		else {
				return 0;
		}
}

IMP::ModelObjectsTemp SurfaceLocalizationRestraint::do_get_inputs() const {
    return plist_;
}


IMPMICOS_END_NAMESPACE
