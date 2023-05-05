#include <math.h>
#include <IMP/core/XYZ.h>
#include <IMP/micos/MembraneInclusionRestraint.h>

IMPMICOS_BEGIN_NAMESPACE

// tries to keep the particles within the membrane thickness (R-r)
MembraneInclusionRestraint::MembraneInclusionRestraint(IMP::ParticlesTemp plist, double R, double r, double sigma) :

		Restraint(plist[0]->get_model(), "MembraneInclusionRestraint %1%"),
		plist_(plist),
		Rsq_(R * R),
		rsq_(r * r),
		sigma_(sigma){}

 /* calculate distance of each particle from origin, which is the center of the base of cylinder */

double MembraneInclusionRestraint::getDistance(IMP::Particle* p) const {
    double x = IMP::core::XYZ(p).get_coordinate(0);  // The coordinates of the particle
    double y = IMP::core::XYZ(p).get_coordinate(1);
    // Calculate the coordinate-wise distance from the center

    double radial = (x * x) + (y * y);  // radial = (euclidean distance from center) ^ 2
    if (radial < rsq_) {  // Continue only if the distance lesser than the inner radius
        double deviation = radial + rsq_ - 2 * sqrt(radial*rsq_);
        // deviation = (sqrt(rsq_) - sqrt(radial))^2
        return fabs(deviation);
    }
    if (radial > Rsq_) {  // Continue only if the distance is more than outer radius
        double deviation = radial + Rsq_ - 2 * sqrt(radial*Rsq_);
        // deviation = (sqrt(radial) - sqrt(Rsq_))^2
        return fabs(deviation);
    }
    else {
        return 0;
    }
}


double MembraneInclusionRestraint::unprotected_evaluate(IMP::DerivativeAccumulator* accum) const {
    double score = 0;
    for (unsigned int i=0; i < plist_.size(); i++){
        score += getDistance(plist_[i]);
        }

    if (accum){};
    return (score/sigma_);
}

IMP::ModelObjectsTemp MembraneInclusionRestraint::do_get_inputs() const {
    return plist_;
}


IMPMICOS_END_NAMESPACE
