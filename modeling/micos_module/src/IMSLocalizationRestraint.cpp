#include <math.h>
#include <IMP/core/XYZ.h>
#include <IMP/micos/IMSLocalizationRestraint.h>

IMPMICOS_BEGIN_NAMESPACE

// Intermembrane space is the lumen (less than r)
IMSLocalizationRestraint::IMSLocalizationRestraint(IMP::ParticlesTemp plist, double r, double sigma) :

		Restraint(plist[0]->get_model(), "IMSLocalizationRestraint %1%"),
		sigma_(sigma),
		rsq_(r * r),
		plist_(plist) {}

 /* calculate distance of each particle from origin, which is the center of the base of cylinder */

double IMSLocalizationRestraint::getDistance(IMP::Particle* p) const {
    double x = IMP::core::XYZ(p).get_coordinate(0);  // The coordinates of the particle
    double y = IMP::core::XYZ(p).get_coordinate(1);
    // Calculate the coordinate-wise distance from the center

    double radial = (x * x) + (y * y);  // radial = (euclidean distance from center) ^ 2
    if (radial > rsq_) {  // Continue only if the distance lesser than the inner radius
        double deviation = radial - rsq_;
        // deviation = squared difference of distance of particle and the inner radius
        return fabs(deviation);
    }
    else {
        return 0;
    }
}


double IMSLocalizationRestraint::unprotected_evaluate(IMP::DerivativeAccumulator* accum) const {
    double score = 0;
    for (unsigned int i=0; i < plist_.size(); i++){
        score += getDistance(plist_[i]);
        }

    if (accum){};
    return (score/sigma_);
}

IMP::ModelObjectsTemp IMSLocalizationRestraint::do_get_inputs() const {
    return plist_;
}


IMPMICOS_END_NAMESPACE
