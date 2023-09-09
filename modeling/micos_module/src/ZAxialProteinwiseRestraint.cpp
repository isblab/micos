#include <cmath>
#include <IMP/core/XYZ.h>
#include <IMP/micos/ZAxialProteinwiseRestraint.h>

IMPMICOS_BEGIN_NAMESPACE

// tries to keep the particles above CJ
ZAxialProteinwiseRestraint::ZAxialProteinwiseRestraint(IMP::ParticlesTemp plist, double lower_bound,double upper_bound, double sigma) :

		Restraint(plist[0]->get_model(), "ZAxialProteinwiseRestraint %1%"),
		plist_(plist),
		lower_bound_(lower_bound),
		upper_bound_(upper_bound),
		sigma_(sigma){}


double ZAxialProteinwiseRestraint::unprotected_evaluate(IMP::DerivativeAccumulator* da) const {
	double score = 0;
	for (unsigned int i=0; i < plist_.size(); i++){
		score += IMP::core::XYZ(plist_[i]).get_coordinate(2);
		        }

	double avg = score/plist_.size();

	if (avg<lower_bound_) {
		double deviation = avg*avg + lower_bound_*lower_bound_ - 2*avg*lower_bound_;
		return fabs(deviation)/sigma_;
			}

	if (avg>upper_bound_) {
		double deviation = avg*avg + upper_bound_*upper_bound_ - 2*avg*upper_bound_;
		return fabs(deviation)/sigma_;
		}

	else{
		return 0;
	}
}

IMP::ModelObjectsTemp ZAxialProteinwiseRestraint::do_get_inputs() const {
    return plist_;
}


IMPMICOS_END_NAMESPACE
