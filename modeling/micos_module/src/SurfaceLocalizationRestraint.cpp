#include <math.h>
#include <IMP/core/XYZ.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/micos/SurfaceLocalizationRestraint.h>
#include "IMP/Model.h"

IMPMICOS_BEGIN_NAMESPACE

SurfaceLocalizationRestraint::SurfaceLocalizationRestraint(IMP::ParticlesTemp plist, double r, double sigma, double allowed_dist) :

		Restraint(plist[0]->get_model(), "SurfaceLocalizationRestraint %1%"),
		plist_(plist),
		rsq_(r * r),
		sigma_(sigma),
		allowed_dist_ ((r - allowed_dist)*(r - allowed_dist)) {}

double SurfaceLocalizationRestraint::getCoordinates(IMP::Particle* p) const {

    if (IMP::core::RigidMember::get_is_setup(p)) {
        IMP::core::RigidMember rb_member(p);
        IMP::core::RigidBody rb = rb_member.get_rigid_body();
        double x = rb.get_coordinates()[0];
        double y = rb.get_coordinates()[1];
        return getDistance(x, y);
        }
        
    else {         
        double x = IMP::core::XYZ(p).get_coordinate(0);  // The coordinates of the particle
	double y = IMP::core::XYZ(p).get_coordinate(1);
	return getDistance(x, y);
	}
   }
	
double SurfaceLocalizationRestraint::getDistance(double x , double y) const{	
	
    double radial = (x * x) + (y * y);
    if (radial < allowed_dist_) {
        double deviation = radial + allowed_dist_ - 2 * sqrt(radial * allowed_dist_);
        return fabs(deviation);
       } 
    else if (radial > rsq_) {
        double deviation = radial + rsq_ - 2 * sqrt(radial * rsq_);
        return fabs(deviation);
        }
    else {
        return 0;
    }
}

double SurfaceLocalizationRestraint::unprotected_evaluate(IMP::DerivativeAccumulator* accum) const {
    double score = 0;
    for (unsigned int i=0; i < plist_.size(); i++){
        score += getCoordinates(plist_[i]);
        }

    if (accum){};
    return (score/sigma_);
}

IMP::ModelObjectsTemp SurfaceLocalizationRestraint::do_get_inputs() const {
    return plist_;
}


IMPMICOS_END_NAMESPACE
