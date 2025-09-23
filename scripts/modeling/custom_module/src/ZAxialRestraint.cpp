#include <math.h>
#include <IMP/core/XYZ.h>
#include <IMP/atom/Mass.h>
#include <IMP/micos/ZAxialRestraint.h>

IMPMICOS_BEGIN_NAMESPACE

// tries to keep the particles above ub
ZAxialRestraint::ZAxialRestraint(IMP::ParticlesTemp plist, double ub, double lb, double sigma, String method) :

		Restraint(plist[0]->get_model(), "ZAxialRestraint %1%"),
		plist_(plist),
		ub_(ub),
		lb_(lb),
		sigma_(sigma),
		method_(method){}


double ZAxialRestraint::getDistance(double z) const {

    if (z < ub_) {  // Continue only if z is above ub centre
        double deviation = z*z + ub_*ub_ -2 *ub_*z;
        return fabs(deviation);
    }
    if (z > lb_) {  // Continue only if z is not outside lb
        double deviation = z*z + lb_*lb_ -2*lb_*z;
        return fabs(deviation);
    }
    else {
        return 0;
    }
}


double ZAxialRestraint::unprotected_evaluate(IMP::DerivativeAccumulator* accum) const {
	double score = 0;
	double total_com = 0;
	double avg = 0;
	double mass =0;
	


	if (method_=="average"){
		for (unsigned int i = 0; i < plist_.size(); i++){
			double com = IMP::atom::Mass(plist_[i]).get_mass();
			double z = IMP::core::XYZ(plist_[i]).get_coordinate(2);
			total_com += com*z;
			mass += com;
		}
	 avg = total_com/mass;
	 score += getDistance(avg);
	}
	
	else if (method_=="domain"){
		std::vector<double> score_list;
		for (unsigned int i = 0; i < plist_.size(); i++){
		double z = IMP::core::XYZ(plist_[i]).get_coordinate(2);
		score_list.push_back(getDistance(z));
		}
		auto it = std::find(score_list.begin(), score_list.end(), 0.0);

    		if (it != score_list.end()) {
        	    score = 0.0;
	       }

	else if (method_=="beadwise"){
		for (unsigned int i = 0; i < plist_.size(); i++) {
			double z = IMP::core::XYZ(plist_[i]).get_coordinate(2);
			double dist = getDistance(z);

			if (dist != 0.0) {  // if any particle is not zero, add the distance from the bounds to the score
				score += dist; 
			}
		}

	       }
	else {
        auto minElement = std::min_element(score_list.begin(), score_list.end());
        if (minElement != score_list.end()) {
            score = *minElement;
            }
    	 }		 
       }	
		else {
			for (unsigned int i = 0; i < plist_.size(); i++){
		double z = IMP::core::XYZ(plist_[i]).get_coordinate(2);
		score += getDistance(z);
	    }
	}

    if (accum){};
    return (score/sigma_);
}

IMP::ModelObjectsTemp ZAxialRestraint::do_get_inputs() const {
    return plist_;
}


IMPMICOS_END_NAMESPACE
