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
		max_limit_ (max_limit * max_limit) {}

double SurfaceLocalizationRestraint::getDistance(IMP::Particle* p) const {

    double x = IMP::core::XYZ(p).get_coordinate(0);  // The coordinates of the particle
    double y = IMP::core::XYZ(p).get_coordinate(1);
    // Calculate the coordinate-wise distance from the center

    double radial = (x * x) + (y * y);  // radial = (euclidean distance from center) ^ 2
    if (radial < max_limit_) {  // Continue only if the distance greater than the max limit
        double deviation = radial + max_limit_ - 2 * sqrt(radial*max_limit_);
        // deviation = (sqrt(radial) - sqrt(rsq_))^2
        return fabs(deviation);
    }
		if (radial > rsq_) {
			double deviation = radial + rsq_ - 2* sqrt(radial*rsq_);

			return fabs(deviation);
		}
    else {
        return 0;
    }
}

double SurfaceLocalizationRestraint::unprotected_evaluate(IMP::DerivativeAccumulator* accum) const {
    double score = 0;
    for (unsigned int i=0; i < plist_.size(); i++){
        score += getDistance(plist_[i]);
        }

    if (accum){};
    return (score/sigma_);
}

IMP::ModelObjectsTemp SurfaceLocalizationRestraint::do_get_inputs() const {
    return plist_;
}


IMPMICOS_END_NAMESPACE
