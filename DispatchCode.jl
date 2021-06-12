using CSV, DataFrames, LinearAlgebra, JuMP, Gurobi

cd("/Users/thodoris/MIT/EnergyProject/RTBM_SCED_Input_Output/")

dat1 = CSV.read("input_resource_energy_offer.txt", DataFrame)
total_id_resources = unique(dat1.ID_RESOURCE)

# The vector total_id_resources contains all the 764 resources


#######################################################################

# Online Flag
dat1 = CSV.read("input_resource_scada.txt", DataFrame)
dat1 = sort!(dat1, :ID_RESOURCE)
OnLineFlag = (dat1.SCADA_MW .!= 0)*1
R = Vector(dat1.ID_RESOURCE)


# Control Status
dat1 = CSV.read("output_resource_commitment_dispatch.txt", DataFrame)
dat1 = sort!(dat1, :ID_RESOURCE)
ControlStatus = Vector(dat1.CMODE)

# EnergyStepPrice, EnergyStepWidth
dat1 = CSV.read("input_resource_energy_offer.txt", DataFrame)
dat1 = sort!(dat1, :ID_RESOURCE)

dat3 = groupby(dat1, :ID_RESOURCE)

width = []
for i in 1:764
    v = []
    push!(v, dat3[i].MW[1])
    for j in 2:length(dat3[i].SEGMENT)
        push!(v, dat3[i].MW[j] - dat3[i].MW[j-1])
    end
    push!(width, v)
end


segments = []
prices = []
for i in 1:764
    push!(segments, Vector(dat3[i].SEGMENT))
    push!(prices, Vector(dat3[i].PRICE))
end

EnergyStepPrice_Segments = segments
EnergyStepPrice_Prices = prices




# Ancillary Service Offers
dat1 = CSV.read("input_resource_ancillary_offer.txt", DataFrame)
dat2 = CSV.read("input_product_requirement.txt", DataFrame)
dat1 = sort!(dat1, :ID_RESOURCE)


# RegUpResOfferPrice
dat3 = filter(row ->row[:PRODUCT] == "REGUP", dat1)
dat4 = filter(row ->row[:ID_PRODUCT] == "REGUP", dat2)
dat3 = dropmissing(dat3)
R_REGUP =  unique(dat3.ID_RESOURCE)
RegUpResOfferPrice = dat3.PRICE .+ dat4.MILEAGE_FACTOR[1]*dat3.MILEAGEPRICE


#RegDnResOfferPrice
dat5 = filter(row -> row[:PRODUCT] == "REGDN", dat1)
dat6 = filter(row -> row[:ID_PRODUCT] == "REGDN", dat2)
dat5 = dropmissing(dat5)
R_REGDN =  unique(dat5.ID_RESOURCE)
RegDnResOfferPrice = dat5.PRICE .+ dat6.MILEAGE_FACTOR[1]*dat5.MILEAGEPRICE

#SpinResOfferPrice, SuppResOfferPrice
dat7 = filter(row -> row[:PRODUCT] == "SPIN", dat1)
dat8 = filter(row -> row[:PRODUCT] == "SUPP", dat1)
dat7 = dat7[:,1:4]
dat8 = dat8[:,1:4]
dat7 = dropmissing(dat7)
dat8 = dropmissing(dat8)
SpinResOfferPrice = dat7.PRICE
SuppResOfferPrice = dat8.PRICE
R_SPIN =  unique(dat7.ID_RESOURCE)
R_SUPP =  unique(dat8.ID_RESOURCE)

# |R_REGUP| = (250)
# |R_REGDN| = (271)
# |R_SPIN| =  (322)
# |R_SUPP| = (288)


# Demand Penalty VRL Curves

# ContResDemandPrice , RegUpResDemandPrice, RegDnResDemandPrice
dat1 = CSV.read("input_demand_penalty_vrl_curves.txt", DataFrame)
dat2 = CSV.read("input_product_requirement.txt", DataFrame)
dat11 = filter(row -> row[:CURVE_TYPE] == "CONTINGENCY_RESERVE", dat1)
dat12 = filter(row -> row[:CURVE_TYPE] == "REGUP", dat1)
dat13 = filter(row -> row[:CURVE_TYPE] == "REGDN", dat1)
market_wide_crr = dat2.REQUIREMENT[1] - dat2.REQUIREMENT[4]
MSSC = market_wide_crr / 1.2
ContResDemandPrice = dat11.PRICE
RegUpResDemandPrice  = dat12.PRICE
RegDnResDemandPrice = dat13.PRICE

ContResDemandStep = [0.2*MSSC, 0.7*MSSC, 1.2*MSSC]



RegUpResDemandStep = [0.1,0.2,0.3,0.5,0.7,1]*dat2.REQUIREMENT[4]
RegDnResDemandStep = [0.1,0.2,0.3,0.5,0.7,1]*dat2.REQUIREMENT[3]

#VRLs

# Penalty Prices
LimitPenPrice = filter(row -> row[:CURVE_TYPE] == "RESOURCE_CAPACITY", dat1).PRICE[1]
RampPenPrice = filter(row -> row[:CURVE_TYPE] == "RESOURCE_RAMP", dat1).PRICE[1]
GlobalPowerBalancePenPrice = filter(row -> row[:CURVE_TYPE] == "GLOBAL_POWER_BALANCE", dat1).PRICE[1]
TransLimitPenPrice = filter(row -> row[:CURVE_TYPE] == "FLOWGATE", dat1).PRICE
SpinResPenPrice = filter(row -> row[:CURVE_TYPE] == "SPINNING_RESERVE", dat1).PRICE[1]


# InitialEnergyOutput
df1 = CSV.read("input_resource_scada.txt", DataFrame)
df1 = sort!(df1, :ID_RESOURCE)
InitialEnergyOutput = Vector(df1.SCADA_MW)


# RegModeFlag
RegModeFlag = (ControlStatus .== 2)*1


# Resource Characteristics

# MaxLimit, MinLimit
dat1 = CSV.read("input_resource_current_properties.txt", DataFrame)
dat1 = sort!(dat1, :ID_RESOURCE)
MaxLimit = Vector(dat1.ECOMAX)
MinLimit = Vector(dat1.ECOMIN)

# if CMODE = 3, then ECOMIN=ECOMAX=Initial Energy Output
# WR.WOLF, Initial 1212.7, ECOMAX = 1200,  RampDownRate = 1MW per minute
# if CMODE = 1 or 2, then if InitialEnergyOutput > MaxLimit then
# InitialEnergyOutput = MaxLimit
for i in 1:764
    if ControlStatus[i] == 3
        MaxLimit[i] = InitialEnergyOutput[i]
        MinLimit[i] = InitialEnergyOutput[i]
    elseif (ControlStatus[i] == 1 || ControlStatus[i] == 2) && (InitialEnergyOutput[i] > MaxLimit[i])
        InitialEnergyOutput[i] = MaxLimit[i]

    end
end

# dict_InitialEnergyOutput["FRONTIER"]
# dict_MaxLimit["FRONTIER"]
# dict_MinLimit["FRONTIER"]

for i in 1:764
    if InitialEnergyOutput[i] < MinLimit[i]
        MinLimit[i] = InitialEnergyOutput[i]
    end
end




# Resource Ancillary Services Qualification

# SpinAvailability, SuppAvailability
dat1 = CSV.read("input_resource_product_qualification.txt", DataFrame)
dat1 = sort!(dat1, :ID_RESOURCE)
SpinAvailability = filter(row -> row[:ID_PRODUCT] == "SPIN", dat1)[:,4]
SuppAvailability = filter(row -> row[:ID_PRODUCT] == "SUPP", dat1)[:,4]


# Resource Energy Ramp Rates
dat1 = CSV.read("output_resource_commitment_dispatch.txt", DataFrame)
dat1 = sort!(dat1, :ID_RESOURCE)
RampUpRate = Vector(dat1.AVG_UP_RR)
RampDownRate = Vector(dat1.AVG_DOWN_RR)
MaxOffLineResponse = min.(RampDownRate*10, MaxLimit)

# Resource Reserve Ramp Rates

# RegRampUPRate, RegRampDownRate, SpinRampRate
dat1 = CSV.read("input_resource_ramp_curve_ancillary.txt", DataFrame)
dat1 = sort!(dat1, :ID_RESOURCE)
df1 = CSV.read("output_resource_product.txt", DataFrame)
df1 = sort!(df1, :ID_RESOURCE)
RegRampUpRate = filter(row -> row[:ID_PRODUCT] == "REGUP", df1).RESOURCE_PRODUCT_RR
RegRampDownRate = filter(row -> row[:ID_PRODUCT] == "REGDN", df1).RESOURCE_PRODUCT_RR
#SpinRampRate = filter(row -> row[:ID_PRODUCT] == "SPIN", df1).RESOURCE_PRODUCT_RR

# change spin ramp rate
df1 = CSV.read("output_resource_commitment_dispatch.txt", DataFrame)
df1 = sort!(df1, :ID_RESOURCE)
SpinRampRate = df1.AVG_UP_RR

dat1 = CSV.read("input_resource_product_qualification.txt", DataFrame)
dat1 = sort!(dat1, :ID_RESOURCE)
dat2 = filter(row -> row[:ID_PRODUCT] == "SPIN", dat1)
SpinResAvailabilityFlag = dat2[:,4]


# Distribution Factor

# DF
dat1 = CSV.read("input_resource_pnode_mapping.txt", DataFrame)
dat2 = CSV.read("input_pnode_network_mapping.txt", DataFrame)
dat1 = sort!(dat1, :ID_RESOURCE)
pnodes = dat1.ID_PNODE

keys1 = []
for p in pnodes
    dat3 = filter(row -> row[:ID_PNODE] == p, dat2)
    push!(keys1, [dat3.ID_SUBSTATION[i]*string(dat3.KV[i])*dat3.ID_NODE[i] for i in 1:size(dat3,1)])
end


aux_list = []
for k in 1:764
    aux_list = vcat(aux_list,keys1[k])
end

I1 = unique(aux_list)
I1 = String.(I1)


DF = zeros(length(I1),length(R))
for (i,r) in enumerate(R)
    for k in keys1[i]
        ind = findfirst(x->x==k, I1)
        DF[ind,i] = 1.0
    end
end


# Load Forecast
# SystemDemandMW
dat1 = CSV.read("input_load_forecast.txt", DataFrame)
SystemDemandMW = sum(dat1.FORECAST)

# Transmission Constraints
# RefBusDF, CAMW
dat1 = CSV.read("input_load_internal.txt", DataFrame)
RefBusDF = (1 / sum(dat1.MW))*dat1.MW

keys2 = []
for i in 1:size(dat1,1)
    push!(keys2, dat1.NAME_SUBSTATION[i]*string(dat1.KV[i])*dat1.ID_ND[i])
end

I2 = keys2
I2 = String.(I2)


# nsi in output_solution_interval_results
dat1 = CSV.read("output_solution_interval_results.txt", DataFrame)
NSI = dat1.NSI[1]


# LossSensitivityFactors
dat1 = CSV.read("input_loss_sensitivities.txt", DataFrame)
dLoss = Vector(dat1.LOSS_SENS)

# LossOffsetFactor
dat1 = CSV.read("input_loss_general_data.txt", DataFrame)
LossOffsetFactor = (dat1.LOSS_OFFSET_MW[1] / dat1.TOTAL_LOSS_MW[1])

# ExternalInjectionMW
# dat1 = CSV.read("input_load_external.txt", DataFrame)
# ExternalInjectionMW = Vector(dat1.MW)
# dat1 = CSV.read("input_generation_external.txt", DataFrame)
# ExternalInjectionMW = Vector(dat1.MW)

#df1 = CSV.read("output_load_external.txt", DataFrame)
#sum(df1.MW)

#df2 = CSV.read("output_generator_external.txt", DataFrame)
#sum(df2.MW)

#sum(df1.MW) - sum(df2.MW)

# Transmission Constraints Shift Factors
dat1 = CSV.read("input_constraint_node_shiftfactors.txt", DataFrame)

R_REGUP_C = [i for i in R if !(i in R_REGUP)]
R_REGDN_C = [i for i in R if !(i in R_REGDN)]
R_SPIN_C = [i for i in R if !(i in R_SPIN)]
R_SUPP_C = [i for i in R if !(i in R_SUPP)]


df1 = CSV.read("output_solution_interval_results.txt", DataFrame)
ModeledLosses = df1.TOTAL_LOSSES

#######################################################################






#######################################################################


# Define The Optimization Model


#ST_1 = [1,2,3,4,5,6]
DI = [1]

m = Model(Gurobi.Optimizer)
set_optimizer_attribute(m, "DualReductions", 0)

# Resource energy costs
dict_EnergyStepPrice_Segments = Dict(i => EnergyStepPrice_Segments[j] for (j,i) in enumerate(R))
dict_EnergyStepPrice_Prices = Dict(i => EnergyStepPrice_Prices[j] for (j,i) in enumerate(R))
@variable(m, ResourceEnergyCosts)
@variable(m, EnergyStepClearing[r in R, di in DI, st in dict_EnergyStepPrice_Segments[r]]>=0)
@variable(m, EnergyDispatchTarget[r in R, di in DI])

#populate EnergyStepPrice
# constraint 1.1.1
@constraint(m, [di in DI], ResourceEnergyCosts == (1/12)*sum(sum(dict_EnergyStepPrice_Prices[r][st]*
            EnergyStepClearing[r, di, st] for st in dict_EnergyStepPrice_Segments[r]) for r in R))

# constraint 1.1.2
@constraint(m, [r in R, di in DI], EnergyDispatchTarget[r,di] == sum(EnergyStepClearing[r,di,st] for st in dict_EnergyStepPrice_Segments[r]))


# Regulation up reserve costs
dict_RegUpResOfferPrice = Dict(i => RegUpResOfferPrice[j] for (j,i) in enumerate(R_REGUP))

@variable(m, RegulationUpReserveCosts)
@variable(m, ClearedRegUpRes[r in R, di in DI]>=0)

# constraint 1.2.1
@constraint(m, [di in DI], RegulationUpReserveCosts == (1/12)*sum(dict_RegUpResOfferPrice[r]*ClearedRegUpRes[r,di] for r in R_REGUP))
@constraint(m, [r in R_REGUP_C, di in DI], ClearedRegUpRes[r,di] == 0.0)

# Regulation down reserve costs
dict_RegDnResOfferPrice = Dict(i => RegDnResOfferPrice[j] for (j,i) in enumerate(R_REGDN))

@variable(m, RegulationDownReserveCosts)
@variable(m, ClearedRegDnRes[r in R, di in DI]>=0)

# constraint 1.3.1
@constraint(m, [di in DI], RegulationDownReserveCosts == (1/12)*sum(dict_RegDnResOfferPrice[r]*ClearedRegDnRes[r,di] for r in R_REGDN))
@constraint(m, [r in R_REGDN_C, di in DI], ClearedRegDnRes[r,di] == 0.0)

# Spinning reserve costs
dict_SpinResOfferPrice = Dict(i => SpinResOfferPrice[j] for (j,i) in enumerate(R_SPIN))

@variable(m, SpinningReserveCosts)
@variable(m, ClearedSpinRes[r in R, di in DI]>=0)

# constraint 1.4.1
@constraint(m, [di in DI], SpinningReserveCosts == (1/12)*sum(dict_SpinResOfferPrice[r]*ClearedSpinRes[r,di] for r in R_SPIN))
@constraint(m, [r in R_SPIN_C, di in DI], ClearedSpinRes[r,di] == 0.0)

# Supplemental reserve costs
dict_SuppResOfferPrice = Dict(i => SuppResOfferPrice[j] for (j,i) in enumerate(R_SUPP))

@variable(m, SupplementalReserveCosts)
@variable(m, ClearedSuppRes[r in R, di in DI]>=0)

# constraint 1.5.1
@constraint(m, [di in DI], SupplementalReserveCosts == (1/12)*sum(dict_SuppResOfferPrice[r]*ClearedSuppRes[r,di] for r in R_SUPP))
@constraint(m, [r in R_SUPP_C, di in DI], ClearedSuppRes[r,di] == 0.0)



# System wide contigency reserve value
ST_CONT = [1,2,3]
#@variable(m, SystemWideContingencyReserveValue)
@variable(m, MWContResStepClearing[st in ST_CONT, di in DI]>=0)

#constraint 1.6.1
#@constraint(m, [di in DI], SystemWideContingencyReserveValue == sum(ContResDemandPrice[st]* MWContResStepClearing[st,di] for st in ST_CONT))

# System wide regulation up reserve value
ST = [1,2,3,4,5,6]
#@variable(m, SystemWideRegulationUpReserveValue)
@variable(m, MWRegUpResStepClearing[st in ST, di in DI]>=0)

#constraint 1.7.1
#@constraint(m, [di in DI], SystemWideRegulationUpReserveValue == sum(RegUpResDemandPrice[st]* MWRegUpResStepClearing[st,di] for st in ST_REG))

# System wide regulation down reserve value
#@variable(m, SystemWideRegulationDownReserveValue)
@variable(m, MWRegDnResStepClearing[st in ST, di in DI]>=0)

#constraint 1.8.1
#@constraint(m, [di in DI], SystemWideRegulationDownReserveValue == sum(RegDnResDemandPrice[st]* MWRegDnResStepClearing[st,di] for st in ST_REG))


# Constraint violation penalty costs

@variable(m, ConstraintViolationPenaltyCosts)
@variable(m, ResourceMaximumLimitViolationCost)
@variable(m, ResourceMinimumLimitViolationCost)
@variable(m, ResourceRampRateViolationCost)
@variable(m, GlobalPowerBalanceViolationCost)
#@variable(m, TransmissionDemandCurveViolationCost)
#
# JuMP.value.(ResourceMaximumLimitViolationCost)
# JuMP.value.(ResourceMinimumLimitViolationCost)
# JuMP.value.(ResourceRampRateViolationCost)
# JuMP.value.(GlobalPowerBalanceViolationCost)

#constraint 1.9.1
@constraint(m, ConstraintViolationPenaltyCosts ==
                 ResourceMaximumLimitViolationCost
               + ResourceMinimumLimitViolationCost
               + ResourceRampRateViolationCost
               + GlobalPowerBalanceViolationCost)
               #+ TransmissionDemandCurveViolationCost)

# Resource maximum limit violation cost
@variable(m, MaxLimitViolation[r in R, di in DI]>=0)

#constraint 1.9.1.1.1
@constraint(m, [di in DI], ResourceMaximumLimitViolationCost == (LimitPenPrice/12)*sum(MaxLimitViolation[r,di] for r in R))

# Resource minimum limit violation cost
@variable(m, MinLimitViolation[r in R, di in DI]>=0)

#constraint 1.9.1.2.1
@constraint(m, [di in DI], ResourceMinimumLimitViolationCost == (LimitPenPrice/12)*sum(MinLimitViolation[r,di] for r in R))


# Resource ramp rate violation cost

@variable(m, EnergyRampUpViolation[r in R, di in DI]>=0)
@variable(m, EnergyRampDownViolation[r in R, di in DI]>=0)
@variable(m, SpinRampViolation[r in R, di in DI]>=0)
@variable(m, RegRampUpViolation[r in R, di in DI]>=0)
@variable(m, RegRampDownViolation[r in R, di in DI]>=0)
#
# maximum([JuMP.value.(EnergyRampDownViolation[r,1]) for r in R])
# maximum([JuMP.value.(EnergyRampUpViolation[r,1]) for r in R])
# maximum([JuMP.value.(SpinRampViolation[r,1]) for r in R])
# maximum([JuMP.value.(RegRampUpViolation[r,1]) for r in R])
# maximum([JuMP.value.(RegRampDownViolation[r,1]) for r in R])

#constraint 1.9.1.3.1
@constraint(m, [di in DI], ResourceRampRateViolationCost == (RampPenPrice/12)*sum(EnergyRampUpViolation[r, di]
        + EnergyRampDownViolation[r, di] + SpinRampViolation[r, di]
        + RegRampUpViolation[r, di] + RegRampDownViolation[r, di]  for r in R))

# Global power balance violation cost
@variable(m, GlobalEnergyShortage[di in DI]>=0)
@variable(m, GlobalEnergySurplus[di in DI]>=0)

#constraint 1.9.4.1.1
@constraint(m, [di in DI], GlobalPowerBalanceViolationCost == (1/12)*GlobalPowerBalancePenPrice*
    (GlobalEnergyShortage[di] + GlobalEnergySurplus[di]))

# Transmission curve violation cost
# Define K
#ST = [1,2,3,4,5]
#K = []
#@variable(m, TransLimitStepViolation[k in K, st in ST, di in DI])

#constraint 1.9.5
#@constraint(m, [di in DI], TransmissionDemandCurveViolationCost == (1/12)*sum(sum(TransLimitPenPrice[st]*TransLimitStepViolation[k, st, di] for k in K)
                           #for st in ST))

# Spinning reserve violation cost
@variable(m, SpinReserveViolation[di in DI]>=0)
@variable(m, SpinningReserveViolationCost)

#constraint 1.9.5
@constraint(m, [di in DI], SpinningReserveViolationCost == (1/12)*SpinResPenPrice * SpinReserveViolation[di])


# Objective
@variable(m, obj_v)

@constraint(m, obj_v ==  ResourceEnergyCosts + RegulationUpReserveCosts + RegulationDownReserveCosts +
                        SpinningReserveCosts + SupplementalReserveCosts +
                        ConstraintViolationPenaltyCosts)

@objective(m, Min, obj_v)

# Constraints

# Energy Step Clearing Constraints

dict_OnLineFlag = Dict(i => OnLineFlag[j] for (j,i) in enumerate(R))
dict_EnergyStepWidth = Dict(i => width[j] for (j,i) in enumerate(R))

#constraint 2.1.1
@constraint(m, [r in R, di in DI, st in dict_EnergyStepPrice_Segments[r]], EnergyStepClearing[r,di,st] <= dict_OnLineFlag[r]*dict_EnergyStepWidth[r][st])

# Energy Dispatch Target Constraints
dict_ControlStatus = Dict(i => ControlStatus[j] for (j,i) in enumerate(R))
dict_InitialEnergyOutput = Dict(i => InitialEnergyOutput[j] for (j,i) in enumerate(R))



#constraint 2.1.2.1, 2.1.2.2, 2.1.2.3
for r in R
    for di in DI
        if (dict_ControlStatus[r]==1.0)||(dict_ControlStatus[r]==2.0)
            @constraint(m, EnergyDispatchTarget[r, di] == sum(EnergyStepClearing[r, di, st] for st in dict_EnergyStepPrice_Segments[r]))
        elseif (dict_ControlStatus[r]==3.0)
            @constraint(m, EnergyDispatchTarget[r, di] == dict_InitialEnergyOutput[r])
        else
            @constraint(m, EnergyDispatchTarget[r, di] == 0.0)
        end
    end
end


# Cleared Regulation Reserve Range Constraints

dict_MinLimit = Dict(i => MinLimit[j] for (j,i) in enumerate(R))
dict_MaxLimit = Dict(i => MaxLimit[j] for (j,i) in enumerate(R))
dict_ControlStatus = Dict(i => ControlStatus[j] for (j,i) in enumerate(R))
dict_RegModeFlag = Dict(i => RegModeFlag[j] for (j,i) in enumerate(R))

#constraint 2.1.3.1
@constraint(m, [r in R, di in DI], ClearedRegUpRes[r, di] + ClearedRegDnRes[r, di] <=
            dict_RegModeFlag[r]*(dict_MaxLimit[r]-dict_MinLimit[r]))



# Cleared spinning reserve constraints



dict_SpinResAvailabilityFlag = Dict(i => SpinResAvailabilityFlag[j] for (j,i) in enumerate(R))

#constraint 2.1.4.1
@constraint(m, [r in R, di in DI], ClearedSpinRes[r, di] <= dict_OnLineFlag[r]*dict_SpinResAvailabilityFlag[r]
                *(dict_MaxLimit[r]-dict_MinLimit[r]))

# Cleared supplemental reserve constraints

dict_MaxOffLineResponse  = Dict(i => MaxOffLineResponse[j] for (j,i) in enumerate(R))
dict_SuppResAvailabilityFlag = Dict(i => SuppAvailability[j] for (j,i) in enumerate(R))

#constraint 2.1.5.1
@constraint(m, [r in R, di in DI], ClearedSpinRes[r,di] <= dict_SuppResAvailabilityFlag[r]*
    (dict_OnLineFlag[r]*(dict_MaxLimit[r]-dict_MinLimit[r])
    + (1-dict_OnLineFlag[r])*dict_MaxOffLineResponse[r]))

# Maximum limit constraints

#constraint 2.1.6.1
@constraint(m, [r in R, di in DI], EnergyDispatchTarget[r,di] + ClearedRegUpRes[r,di] +
        ClearedSpinRes[r,di] + ClearedSuppRes[r,di] <=  dict_OnLineFlag[r]*dict_MaxLimit[r] +
        (1-dict_OnLineFlag[r])*dict_MaxOffLineResponse[r] + MaxLimitViolation[r,di])

# Minimum limit constraints

#constraint 2.1.7.1
@constraint(m, [r in R, di in DI], EnergyDispatchTarget[r,di] - ClearedRegDnRes[r,di]
    >= dict_OnLineFlag[r]*dict_MinLimit[r] - MinLimitViolation[r, di])


# Energy ramp up constraints
dict_RampUpRate = Dict(i => RampUpRate[j] for (j,i) in enumerate(R))
dict_RampDownRate = Dict(i => RampDownRate[j] for (j,i) in enumerate(R))
dict_RegRampUpRate = Dict(i => RegRampUpRate[j] for (j,i) in enumerate(R))
dict_RegRampDownRate = Dict(i => RegRampDownRate[j] for (j,i) in enumerate(R))

# R[argmax([JuMP.value.(EnergyRampUpViolation[r,1]) for r in R])]
#
# dict_InitialEnergyOutput["FRONTIER"]
# JuMP.value.(EnergyDispatchTarget["FRONTIER",1])
# dict_RegRampUpRate["FRONTIER"]
# dict_RampUpRate["FRONTIER"]
# dict_RegModeFlag["FRONTIER"]
#
# dict_ControlStatus["FRONTIER"]

#constraint 2.1.8.1
@constraint(m, [r in R, di in DI], EnergyDispatchTarget[r,di] - dict_InitialEnergyOutput[r] +
              ClearedRegUpRes[r, di] <= 5*dict_RegRampUpRate[r] +
              5*(1-dict_RegModeFlag[r])*(dict_RampUpRate[r]-dict_RegRampUpRate[r])
              + EnergyRampUpViolation[r, di])


# Energy ramp down constraints
#constraint 2.1.9.1
@constraint(m, [r in R, di in DI], EnergyDispatchTarget[r,di] - dict_InitialEnergyOutput[r] +
               ClearedRegDnRes[r,di] >= (-5)*dict_RegRampDownRate[r] -
               5*(1-dict_RegModeFlag[r])*(dict_RampDownRate[r]-dict_RegRampDownRate[r])
               + EnergyRampDownViolation[r, di])

# Spinning ramp constraints


dict_SpinRampRate = Dict(i => SpinRampRate[j] for (j,i) in enumerate(R))
ContResDeployTime = 10

#constraint 2.1.10.1
@constraint(m, [r in R, di in DI], ClearedSpinRes[r,di] <= ContResDeployTime
              * (dict_SpinRampRate[r]-(1/5)*(ClearedRegDnRes[r,di]
               + EnergyDispatchTarget[r,di]-dict_InitialEnergyOutput[r]))
               + SpinRampViolation[r,di])


# Regulation ramp up constraints
#constraint 2.1.11.1
@constraint(m, [r in R, di in DI], ClearedRegUpRes[r,di] <= 5*dict_RegRampUpRate[r] + RegRampUpViolation[r,di])


# Regulation ramp down constraints
#constraint 2.1.12.1
@constraint(m, [r in R, di in DI], ClearedRegDnRes[r,di]<= 5*dict_RegRampDownRate[r] + RegRampDownViolation[r,di])


# extra
######### Global Power Balance Constraints #########


# Bus energy injection constraints
# @variable(m, PI[i in 1:length(I1), d in DI])
#
# #constraint 2.2.1.1
# @constraint(m, [i in 1:length(I1), di in DI], PI[i, di] == sum(DF[i, j]*EnergyDispatchTarget[r, di] for (j,r) in enumerate(R)))
#
#
# # Bus energy withdrawal constraints
# @variable(m, PW[i in 1:length(I2), di in DI])
#
# #constraint 2.2.2.1
# @constraint(m, [i in 1:length(I2), di in DI], PW[i, di] == RefBusDF[i]* SystemDemandMW)
#


# Modeled losses constraints # LossOffesetMW variable or input
 #@variable(m, ModeledLosses[di in DI])  

#constraint 2.2.3.1
 #@constraint(m, [di in DI], ModeledLosses[di] == sum(dLoss[i, di]*(PI[i, di]-PW[i, di]+ExternalInjection[i, di]) for i in I) + LossOffsetFactor)

# constraint 2.2.4.1
# @constraint(m, [di in DI], sum(PI[i, di] for i in 1:length(I1))- NSI + GlobalEnergyShortage[di] 
#                     == sum(PW[i, di] for i in 1:length(I2)) + ModeledLosses[1] + GlobalEnergySurplus[di])


#extra
# constraint 2.2.4.1
 @constraint(m, [di in DI], sum(EnergyDispatchTarget[r,di] for r in R)- NSI + GlobalEnergyShortage[di] 
                    == SystemDemandMW + GlobalEnergySurplus[di])



#constraint 2.2.5.3
@constraint(m, [di in DI, st in ST], MWRegUpResStepClearing[st,di] <= RegUpResDemandStep[st])


#constraint 2.2.5.4
@constraint(m, [di in DI, st in ST], MWRegDnResStepClearing[st,di] <= RegDnResDemandStep[st])

#constraint 2.2.6.1 
#@constraint(m, [di in DI], sum((ClearedRegUpRes[r, di]+ClearedSpinRes[r, di]) for r in R) >= SpinRegReq[di] + SpinReserveViolation[di])

#constraint 2.2.7.1 
@constraint(m, [di in DI], sum((ClearedRegUpRes[r, di]+ClearedSpinRes[r, di]+ClearedSuppRes[r, di]) for r in R)  >=
                              sum(MWRegUpResStepClearing[st, di] for st in ST)+ sum(MWContResStepClearing[st,di] for st in ST_CONT))

#constraint 2.2.7.2
@constraint(m, [di in DI, st in ST_CONT], MWContResStepClearing[st,di] <= ContResDemandStep[st])


optimize!(m)


JuMP.value.(ResourceEnergyCosts)

JuMP.value.(RegulationUpReserveCosts)

JuMP.value.(RegulationDownReserveCosts)

JuMP.value.(SpinningReserveCosts)

JuMP.value.(SupplementalReserveCosts)

JuMP.value.(ConstraintViolationPenaltyCosts)



# Comparisons with Output

df1 = CSV.read("output_resource_commitment_dispatch.txt", DataFrame)
df2 = sort!(df1, :ID_RESOURCE)
df5 = hcat(JuMP.value.(EnergyDispatchTarget),df2.DISPATCH_MW)
entries_dif = [i for i in 1:size(df5)[1] if abs(df5[i,1] - df5[i,2]) >= 1]
sum(abs.(df5[:,1] .- df5[:,2]))

EnergyDispatchTarget_Final = df5


df1 = CSV.read("output_resource_product.txt", DataFrame)
df2 = filter(row -> row[:ID_PRODUCT] == "REGUP", df1)
df3 = filter(row -> row[:ID_PRODUCT] == "REGDN", df1)
df4 = filter(row -> row[:ID_PRODUCT] == "SPIN", df1)
df5 = filter(row -> row[:ID_PRODUCT] == "SUPP", df1)
ClearedRegUp_Final = hcat(JuMP.value.(ClearedRegUpRes),df2.PRODUCT_CLEARED)
ClearedRegDn_Final = hcat(JuMP.value.(ClearedRegUpRes),df3.PRODUCT_CLEARED)
ClearedSpinRes_Final = hcat(JuMP.value.(ClearedSpinRes),df4.PRODUCT_CLEARED)
ClearedSuppRes_Final = hcat(JuMP.value.(ClearedSuppRes),df5.PRODUCT_CLEARED)

df6 = EnergyDispatchTarget_Final
entries_dif = [i for i in 1:size(df6)[1] if abs(df6[i,1] - df6[i,2]) >= 0.001]
sum([abs(df6[i,1] - df6[i,2]) for i in 1:size(df6)[1]]) / size(df6)[1]
