################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../AdvancedEstimator.cpp \
../Alligner.cpp \
../ChromosomeInitializer.cpp \
../ClusterComponent.cpp \
../ClusterEstimator.cpp \
../Controller.cpp \
../DirichletDistribution.cpp \
../DirichletMixture.cpp \
../Downloader.cpp \
../Estimator.cpp \
../Functions.cpp \
../Operators.cpp \
../ParameterHandler.cpp \
../RNAread.cpp \
../Run.cpp \
../SimpleEstimator.cpp \
../Simulator.cpp \
../main.cpp \
../main_Test.cpp \
../myRandomEngine.cpp 

O_SRCS += \
../AdvancedEstimator.o \
../Alligner.o \
../ChromosomeInitializer.o \
../ClusterComponent.o \
../ClusterEstimator.o \
../Controller.o \
../DirichletDistribution.o \
../DirichletMixture.o \
../Downloader.o \
../Estimator.o \
../Functions.o \
../Operators.o \
../ParameterHandler.o \
../RNAread.o \
../Run.o \
../SimpleEstimator.o \
../Simulator.o \
../main.o \
../myRandomEngine.o 

OBJS += \
./AdvancedEstimator.o \
./Alligner.o \
./ChromosomeInitializer.o \
./ClusterComponent.o \
./ClusterEstimator.o \
./Controller.o \
./DirichletDistribution.o \
./DirichletMixture.o \
./Downloader.o \
./Estimator.o \
./Functions.o \
./Operators.o \
./ParameterHandler.o \
./RNAread.o \
./Run.o \
./SimpleEstimator.o \
./Simulator.o \
./main.o \
./main_Test.o \
./myRandomEngine.o 

CPP_DEPS += \
./AdvancedEstimator.d \
./Alligner.d \
./ChromosomeInitializer.d \
./ClusterComponent.d \
./ClusterEstimator.d \
./Controller.d \
./DirichletDistribution.d \
./DirichletMixture.d \
./Downloader.d \
./Estimator.d \
./Functions.d \
./Operators.d \
./ParameterHandler.d \
./RNAread.d \
./Run.d \
./SimpleEstimator.d \
./Simulator.d \
./main.d \
./main_Test.d \
./myRandomEngine.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


