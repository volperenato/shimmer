//-------------------------------------------------------------------------------------------------------
//  Shimmer.h
//  BranchReverb plugin based on the Freeverb scheme
//  Customized by Renato Volpe on 07/02/2022.
//
//-------------------------------------------------------------------------------------------------------

#pragma once
#include <stdio.h>
#include "FDN.h"
#include "PitchShifter.h"
#include "vst-2.4-sdk/vstsdk2.4/public.sdk/source/vst2.x/audioeffectx.h"
#include <math.h>
#include "PSMVocoder.h"

using namespace std;

// declare enum for reverb's parameters
enum EfxParameter {
	Param_mix = 0,
	Param_roomSize,
	Param_decay,
	Param_shimmer,
	Param_intervals,
	Param_damping,
	Param_space,
	Param_modDepth,
	Param_modRate,
	Param_lpf,
	Param_hpf,
	Param_Count
};

struct ShimmerParameters
{
	float mix, size, decay, shimmer, interval, damping, rate, depth, space, lpf, hpf;
};

// Declare class BranchReverb
class Shimmer;

// Declare class to handle reverb's presets
class ShimmerPresets {
	friend class Shimmer;
private:
	// Shimmer User Parameters
	float shim_mix, shim_roomSize, shim_shimmer, shim_intervals, shim_decay, shim_damping, shim_space, shim_modRate, shim_modDepth, shim_lpf, shim_hpf; 
	char name[24];
};

//-------------------------------------------------------------------------------------------------------
class Shimmer : public AudioEffectX {	

	// Initialize ShimmerPresets instance
	ShimmerPresets* shim_presets;
	
	// Shimmer User Parameters
	float shim_mix, shim_roomSize, shim_shimmer, shim_intervals, shim_decay, shim_damping, shim_space, shim_modRate, shim_modDepth, shim_lpf, shim_hpf;

	// FDN reverb
	FDN* BranchReverb;
	FDN* MasterReverb;

	// Pitch Shifters
	PSMVocoder* PitchShift_1oct;
	PSMVocoder* PitchShift_2oct;

	// Internal quantities
	float _wet, _dry, _mixP1, _mixP2;	

private:

	void InitPlugin();
	void InitPresets();
	void updateMix();
	void updateMixPitchShifters(float pitch2);

public:

	Shimmer(audioMasterCallback audioMaster);
	~Shimmer();

	// Processing
	virtual void processReplacing(float** inputs, float** outputs, VstInt32 sampleFrames) override;
	virtual float getParameter(VstInt32 index) override;
	virtual void setParameter(VstInt32 index, float value) override;
	virtual bool getEffectName(char* name) override;
	virtual bool getVendorString(char* name) override;
	virtual void getParameterLabel(VstInt32 index, char* label) override;
	virtual void getParameterDisplay(VstInt32 index, char* text) override;
	virtual void getParameterName(VstInt32 index, char* text) override;
	virtual void setSampleRate(float sampleRate) override;
	virtual void setProgram(VstInt32 program) override;
	virtual void getProgramName(char* name) override;
	virtual bool getProgramNameIndexed(VstInt32 category, VstInt32 index, char* text) override;
	virtual VstInt32 getChunk(void** data, bool isPreset) override;
	virtual VstInt32 setChunk(void* data, VstInt32 byteSize, bool isPreset) override;
};





