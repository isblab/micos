#include <math.h>
#include <IMP/core/XYZ.h>
#include <IMP/micos/SurfaceLocalizationRestraint.h>

IMPMICOS_BEGIN_NAMESPACE

// Intermembrane space is the lumen (less than r)
SurfaceLocalizationRestraint::SurfaceLocalizationRestraint(IMP::, double r, double sigma) :

		Restraint(plist[0]->get_model(), "SurfaceLocalizationRestraint %1%"),
		sigma_(sigma),
		rsq_(r * r),
		plist_(plist) {}

 /* calculate distance of each particle from origin, which is the center of the base of cylinder */

double SurfaceLocalizationRestraint::getDistance(IMP::) const {
    double x = IMP::core::RigidBody.get_coordinate(0);  // The global coordinates of the rb
    double y = IMP::core::RigidBody.get_coordinate(1);
    // Calculate the coordinate-wise distance from the center

    double radial = (x * x) + (y * y);  // radial = (euclidean distance from center) ^ 2
    if (sqrt(rsq_) - sqrt(radial)) < 5 {  // Continue only if the distance lesser than the 5A, an arbitrary value
        double deviation = rsq_ -radial_;
        // deviation = squared difference of distance of particle and the inner radius
        return fabs(deviation);
    }
    else {
        return 0;
    }
}


double SurfaceLocalizationRestraint::unprotected_evaluate(IMP::DerivativeAccumulator* accum) const {
    double score = getDistance(rb_);

    if (accum){};
    return (score/sigma_);
}

IMP::ModelObjectsTemp SurfaceLocalizationRestraint::do_get_inputs() const {
    return rb_;
}


IMPMICOS_END_NAMESPACE
