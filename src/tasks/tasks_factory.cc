#include "tasks_factory.h"

using namespace Core;


namespace Tasks
{


PtrTask TaskFactory::create(FILTER_TYPE ftype, TASK_ID id, APPROX_TYPE atype)
{
    switch (ftype) {
    case FILTER_TYPE::Continuous:
        return createContinuous(id, atype);
    case FILTER_TYPE::ContinuousDiscrete:
        return createContinuousDiscrete(id, atype);
    case FILTER_TYPE::Discrete:
        return createDiscrete(id, atype);
    }
    return PtrTask(nullptr);
}

PtrTask TaskFactory::createContinuous(TASK_ID id, APPROX_TYPE type)
{
    switch (type) {
    case APPROX_TYPE::Linear:
        switch (id) {
        case TASK_ID::Landing:
            return PtrTask(new Tasks::Continuous::LandingLinear);
        case TASK_ID::VanDerPol:
            return PtrTask(new Tasks::Continuous::VanDerPolLinear);
        }
    case APPROX_TYPE::Gauss:
        switch (id) {
        case TASK_ID::Landing:
            return PtrTask(new Tasks::Continuous::LandingGauss);
        case TASK_ID::VanDerPol:
            return PtrTask(new Tasks::Continuous::VanDerPolGauss);
        }
    }
    return PtrTask(nullptr);
}

PtrTask TaskFactory::createContinuousDiscrete(TASK_ID id, APPROX_TYPE type)
{
    switch (type) {
    case APPROX_TYPE::Linear:
        switch (id) {
        case TASK_ID::Landing:
            return PtrTask(new Tasks::ContinuousDiscrete::LandingLinear);
        case TASK_ID::VanDerPol:
            return PtrTask(new Tasks::ContinuousDiscrete::VanDerPolLinear);
        }
    case APPROX_TYPE::Gauss:
        switch (id) {
        case TASK_ID::Landing:
            return PtrTask(new Tasks::ContinuousDiscrete::LandingGauss);
        case TASK_ID::VanDerPol:
            return PtrTask(new Tasks::ContinuousDiscrete::VanDerPolGauss);
        }
    }
    return PtrTask(nullptr);
}

PtrTask TaskFactory::createDiscrete(TASK_ID id, APPROX_TYPE type)
{
    switch (type) {
    case APPROX_TYPE::Linear:
        switch (id) {
        case TASK_ID::Landing:
            return PtrTask(new Tasks::Discrete::LandingLinear);
        case TASK_ID::VanDerPol:
            return PtrTask(nullptr); // WARNING
        }
    case APPROX_TYPE::Gauss:
        switch (id) {
        case TASK_ID::Landing:
            return PtrTask(new Tasks::Discrete::LandingGauss);
        case TASK_ID::VanDerPol:
            return PtrTask(nullptr); // WARNING
        }
    }
    return PtrTask(nullptr);
}


} // end Tasks
