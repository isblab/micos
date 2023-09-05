#include <math.h>
#include <IMP/core/XYZ.h>
#include <IMP/micos/ZAxialRestraint.h>

IMPMICOS_BEGIN_NAMESPACE

// tries to keep the particles above CJ
ZAxialRestraint::ZAxialRestraint(IMP::ParticlesTemp plist, double cj, double om, double sigma) :

		Restraint(plist[0]->get_model(), "ZAxialRestraint %1%"),
		plist_(plist),
		cj_(cj),
		om_(om),
		sigma_(sigma){}


double ZAxialRestraint::getDistance(IMP::Particle* p) const {
    double z = IMP::core::XYZ(p).get_coordinate(2);  // The coordinates of the particle


    if (z > cj_) {  // Continue only if z is above cj centre
        double deviation = z*z + cj_*cj_ -2 *cj_*z;
        return fabs(deviation);
    }
    if (z < om_) {  // Continue only if z is not outside OM
        double deviation = z*z + om_*om_ -2*om_*z;
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
