#include <math.h>
#include <IMP/core/XYZ.h>
#include <IMP/micos/ZAxialRestraint.h>

IMPMICOS_BEGIN_NAMESPACE

// tries to keep the particles above CJ
ZAxialRestraint::ZAxialRestraint(IMP::ParticlesTemp plist, double sigma) :

		Restraint(plist[0]->get_model(), "ZAxialRestraint %1%"),
		plist_(plist),
		sigma_(sigma){}


double ZAxialRestraint::getDistance(IMP::Particle* p) const {
    double z = IMP::core::XYZ(p).get_coordinate(2);  // The coordinates of the particle


    if (z > 0) {  // Continue only if z is above cj centre
        double deviation = z-0;
        return fabs(deviation);
    }
    if (z < -50) {  // Continue only if z is not outside OM
        double deviation = -50-z;
        return fabs(deviation);
    }
    else {
        return 0;
    }
}


double ZAxialRestraint::unprotected_evaluate(IMP::DerivativeAccumulator* accum) const {
    double score = 0;
    for (unsigned int i=0; i < plist_.size(); i++){
        score += getDistance(plist_[i]);
        }

    if (accum){};
    return (score/sigma_);
}

IMP::ModelObjectsTemp ZAxialRestraint::do_get_inputs() const {
    return plist_;
}


IMPMICOS_END_NAMESPACE
