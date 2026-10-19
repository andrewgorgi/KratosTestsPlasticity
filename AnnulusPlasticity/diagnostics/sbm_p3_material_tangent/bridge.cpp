#include <pybind11/pybind11.h>
#include "includes/define_python.h"
#include "custom_constitutive/small_strains/plasticity/small_strain_j2_plasticity_plane_strain_2d_three_components.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
namespace Kratos {
class DiagnosticNumericalJ2 : public SmallStrainJ2PlasticityPlaneStrain2DThreeComponents {
public:
    KRATOS_CLASS_POINTER_DEFINITION(DiagnosticNumericalJ2);
    using Base=SmallStrainJ2PlasticityPlaneStrain2DThreeComponents;
    DiagnosticNumericalJ2(int Order,double H):mOrder(Order),mH(H){}
    DiagnosticNumericalJ2(const DiagnosticNumericalJ2& Other):Base(Other),mOrder(Other.mOrder),mH(Other.mH){}
    ConstitutiveLaw::Pointer Clone() const override {return Kratos::make_shared<DiagnosticNumericalJ2>(*this);}
protected:
    void CalculateStressResponse(Parameters& rValues,Vector& rPlastic,double& rAlpha) override {
        const bool stress_flag=rValues.GetOptions().Is(COMPUTE_STRESS);
        if(rValues.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR))rValues.GetOptions().Set(COMPUTE_STRESS,true);
        Base::CalculateStressResponse(rValues,rPlastic,rAlpha);
        if(rValues.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR)) {
            Parameters values=rValues;
            Properties properties=rValues.GetMaterialProperties();
            if(mH>0.)properties.SetValue(PERTURBATION_SIZE,mH);
            values.SetMaterialProperties(properties);
            TangentOperatorCalculatorUtility::CalculateTangentTensorSmallDeformationProvidedStrain(values,this,StressMeasure_Cauchy,true,mOrder);
        }
        rValues.GetOptions().Set(COMPUTE_STRESS,stress_flag);
    }
private:int mOrder;double mH;
};
}
PYBIND11_MODULE(sbm_tangent_bridge,m) {
    namespace py=pybind11;using namespace Kratos;
    py::class_<DiagnosticNumericalJ2,DiagnosticNumericalJ2::Pointer,ConstitutiveLaw>(m,"NumericalJ2").def(py::init<int,double>());
    m.def("numerical",[](ConstitutiveLaw::Parameters& original,ConstitutiveLaw& law,int order,double h){
        ConstitutiveLaw::Parameters values=original;
        Properties properties=original.GetMaterialProperties();
        if(h>0.)properties.SetValue(PERTURBATION_SIZE,h);
        values.SetMaterialProperties(properties);
        TangentOperatorCalculatorUtility::CalculateTangentTensorSmallDeformationProvidedStrain(values,&law,ConstitutiveLaw::StressMeasure_Cauchy,true,order);
        return Matrix(values.GetConstitutiveMatrix());
    });
}
