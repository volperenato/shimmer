; Listing generated by Microsoft (R) Optimizing Compiler Version 19.29.30145.0 

include listing.inc

INCLUDELIB OLDNAMES

PUBLIC	??_R1A@?0A@EJ@Delay@@8				; Delay::`RTTI Base Class Descriptor at (0,-1,0,73)'
PUBLIC	??_R1A@?0A@EA@ModDelay@@8			; ModDelay::`RTTI Base Class Descriptor at (0,-1,0,64)'
PUBLIC	??_R0?AVModDelay@@@8				; ModDelay `RTTI Type Descriptor'
PUBLIC	??_R3ModDelay@@8				; ModDelay::`RTTI Class Hierarchy Descriptor'
PUBLIC	??_R1A@?0A@EN@CombFilter@@8			; CombFilter::`RTTI Base Class Descriptor at (0,-1,0,77)'
PUBLIC	??_7ModDelay@@6B@				; ModDelay::`vftable'
PUBLIC	??_R2ModDelay@@8				; ModDelay::`RTTI Base Class Array'
PUBLIC	??_R4ModDelay@@6B@				; ModDelay::`RTTI Complete Object Locator'
;	COMDAT ??_R4ModDelay@@6B@
rdata$r	SEGMENT
??_R4ModDelay@@6B@ DD 01H				; ModDelay::`RTTI Complete Object Locator'
	DD	00H
	DD	00H
	DD	imagerel ??_R0?AVModDelay@@@8
	DD	imagerel ??_R3ModDelay@@8
	DD	imagerel ??_R4ModDelay@@6B@
rdata$r	ENDS
;	COMDAT ??_R2ModDelay@@8
rdata$r	SEGMENT
??_R2ModDelay@@8 DD imagerel ??_R1A@?0A@EA@ModDelay@@8	; ModDelay::`RTTI Base Class Array'
	DD	imagerel ??_R1A@?0A@EN@CombFilter@@8
	DD	imagerel ??_R1A@?0A@EJ@Delay@@8
	ORG $+3
rdata$r	ENDS
;	COMDAT ??_7ModDelay@@6B@
CONST	SEGMENT
??_7ModDelay@@6B@ DQ FLAT:??_R4ModDelay@@6B@		; ModDelay::`vftable'
	DQ	FLAT:?init@Delay@@UEAAXMH@Z
	DQ	FLAT:?setSampleRate@ModDelay@@UEAAXH@Z
	DQ	FLAT:?processAudio@ModDelay@@UEAAMM@Z
	DQ	FLAT:?init@ModDelay@@UEAAXHW4ModulationType@@MM@Z
CONST	ENDS
;	COMDAT ??_R1A@?0A@EN@CombFilter@@8
rdata$r	SEGMENT
??_R1A@?0A@EN@CombFilter@@8 DD imagerel ??_R0?AVCombFilter@@@8 ; CombFilter::`RTTI Base Class Descriptor at (0,-1,0,77)'
	DD	01H
	DD	00H
	DD	0ffffffffH
	DD	00H
	DD	04dH
	DD	imagerel ??_R3CombFilter@@8
rdata$r	ENDS
;	COMDAT ??_R3ModDelay@@8
rdata$r	SEGMENT
??_R3ModDelay@@8 DD 00H					; ModDelay::`RTTI Class Hierarchy Descriptor'
	DD	00H
	DD	03H
	DD	imagerel ??_R2ModDelay@@8
rdata$r	ENDS
;	COMDAT ??_R0?AVModDelay@@@8
data$rs	SEGMENT
??_R0?AVModDelay@@@8 DQ FLAT:??_7type_info@@6B@		; ModDelay `RTTI Type Descriptor'
	DQ	0000000000000000H
	DB	'.?AVModDelay@@', 00H
data$rs	ENDS
;	COMDAT ??_R1A@?0A@EA@ModDelay@@8
rdata$r	SEGMENT
??_R1A@?0A@EA@ModDelay@@8 DD imagerel ??_R0?AVModDelay@@@8 ; ModDelay::`RTTI Base Class Descriptor at (0,-1,0,64)'
	DD	02H
	DD	00H
	DD	0ffffffffH
	DD	00H
	DD	040H
	DD	imagerel ??_R3ModDelay@@8
rdata$r	ENDS
;	COMDAT ??_R1A@?0A@EJ@Delay@@8
rdata$r	SEGMENT
??_R1A@?0A@EJ@Delay@@8 DD imagerel ??_R0?AVDelay@@@8	; Delay::`RTTI Base Class Descriptor at (0,-1,0,73)'
	DD	00H
	DD	00H
	DD	0ffffffffH
	DD	00H
	DD	049H
	DD	imagerel ??_R3Delay@@8
?_Fake_alloc@std@@3U_Fake_allocator@1@B	ORG $+1		; std::_Fake_alloc
PUBLIC	?processAudio@ModDelay@@UEAAMM@Z		; ModDelay::processAudio
PUBLIC	?setSampleRate@ModDelay@@UEAAXH@Z		; ModDelay::setSampleRate
PUBLIC	?setModType@ModDelay@@QEAAXW4ModulationType@@@Z	; ModDelay::setModType
PUBLIC	?setModDepth@ModDelay@@QEAAXM@Z			; ModDelay::setModDepth
PUBLIC	?init@ModDelay@@UEAAXHW4ModulationType@@MM@Z	; ModDelay::init
PUBLIC	??_GLFO@@QEAAPEAXI@Z				; LFO::`scalar deleting destructor'
PUBLIC	??1ModDelay@@QEAA@XZ				; ModDelay::~ModDelay
;	COMDAT pdata
pdata	SEGMENT
$pdata$?processAudio@ModDelay@@UEAAMM@Z DD imagerel $LN7
	DD	imagerel $LN7+113
	DD	imagerel $unwind$?processAudio@ModDelay@@UEAAMM@Z
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$?setSampleRate@ModDelay@@UEAAXH@Z DD imagerel $LN4
	DD	imagerel $LN4+45
	DD	imagerel $unwind$?setSampleRate@ModDelay@@UEAAXH@Z
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$?init@ModDelay@@UEAAXHW4ModulationType@@MM@Z DD imagerel $LN16
	DD	imagerel $LN16+320
	DD	imagerel $unwind$?init@ModDelay@@UEAAXHW4ModulationType@@MM@Z
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$??_GLFO@@QEAAPEAXI@Z DD imagerel $LN6
	DD	imagerel $LN6+52
	DD	imagerel $unwind$??_GLFO@@QEAAPEAXI@Z
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$??1ModDelay@@QEAA@XZ DD imagerel $LN13
	DD	imagerel $LN13+81
	DD	imagerel $unwind$??1ModDelay@@QEAA@XZ
;	COMDAT xdata
xdata	SEGMENT
$unwind$??1ModDelay@@QEAA@XZ DD 040a01H
	DD	06340aH
	DD	07006320aH
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$??_GLFO@@QEAAPEAXI@Z DD 040a01H
	DD	06340aH
	DD	07006320aH
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$?init@ModDelay@@UEAAXHW4ModulationType@@MM@Z DD 061801H
	DD	026818H
	DD	08340aH
	DD	07006520aH
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$?setSampleRate@ModDelay@@UEAAXH@Z DD 040a01H
	DD	06340aH
	DD	07006320aH
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$?processAudio@ModDelay@@UEAAMM@Z DD 040e01H
	DD	02680eH
	DD	030025206H
; Function compile flags: /Ogspy
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\CombFilter.cpp
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\CombFilter.cpp
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp
;	COMDAT ??1ModDelay@@QEAA@XZ
_TEXT	SEGMENT
this$ = 48
??1ModDelay@@QEAA@XZ PROC				; ModDelay::~ModDelay, COMDAT

; 14   : {

$LN13:
	mov	QWORD PTR [rsp+8], rbx
	push	rdi
	sub	rsp, 32					; 00000020H

; 15   : 	delete mdly_lfo;

	mov	rdi, QWORD PTR [rcx+72]
	lea	rax, OFFSET FLAT:??_7ModDelay@@6B@
	mov	QWORD PTR [rcx], rax
	mov	rbx, rcx
	test	rdi, rdi
	je	SHORT $LN6@ModDelay
	mov	rcx, rdi
	call	??1LFO@@QEAA@XZ				; LFO::~LFO
	mov	edx, 72					; 00000048H
	mov	rcx, rdi
	call	??3@YAXPEAX_K@Z				; operator delete
$LN6@ModDelay:
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\CombFilter.cpp

; 14   : CombFilter::~CombFilter() {}

	lea	rax, OFFSET FLAT:??_7CombFilter@@6B@
	mov	rcx, rbx
	mov	QWORD PTR [rbx], rax
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp

; 16   : }

	mov	rbx, QWORD PTR [rsp+48]
	add	rsp, 32					; 00000020H
	pop	rdi
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\CombFilter.cpp

; 14   : CombFilter::~CombFilter() {}

	jmp	??1Delay@@QEAA@XZ			; Delay::~Delay
??1ModDelay@@QEAA@XZ ENDP				; ModDelay::~ModDelay
_TEXT	ENDS
; Function compile flags: /Ogspy
;	COMDAT ??_GLFO@@QEAAPEAXI@Z
_TEXT	SEGMENT
this$ = 48
__flags$ = 56
??_GLFO@@QEAAPEAXI@Z PROC				; LFO::`scalar deleting destructor', COMDAT
$LN6:
	mov	QWORD PTR [rsp+8], rbx
	push	rdi
	sub	rsp, 32					; 00000020H
	mov	ebx, edx
	mov	rdi, rcx
	call	??1LFO@@QEAA@XZ				; LFO::~LFO
	test	bl, 1
	je	SHORT $LN2@scalar
	mov	edx, 72					; 00000048H
	mov	rcx, rdi
	call	??3@YAXPEAX_K@Z				; operator delete
$LN2@scalar:
	mov	rbx, QWORD PTR [rsp+48]
	mov	rax, rdi
	add	rsp, 32					; 00000020H
	pop	rdi
	ret	0
??_GLFO@@QEAAPEAXI@Z ENDP				; LFO::`scalar deleting destructor'
_TEXT	ENDS
; Function compile flags: /Ogspy
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\Delay.cpp
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp
;	COMDAT ?init@ModDelay@@UEAAXHW4ModulationType@@MM@Z
_TEXT	SEGMENT
this$ = 64
sampleRate$ = 72
modType$ = 80
modRate$ = 88
modDepth$ = 96
?init@ModDelay@@UEAAXHW4ModulationType@@MM@Z PROC	; ModDelay::init, COMDAT

; 19   : {

$LN16:
	mov	QWORD PTR [rsp+8], rbx
	push	rdi
	sub	rsp, 48					; 00000030H
	mov	rbx, rcx

; 49   : 	mdly_modType = modType;

	mov	DWORD PTR [rcx+96], r8d

; 50   : 	switch (mdly_modType)

	xor	ecx, ecx
	movaps	XMMWORD PTR [rsp+32], xmm6

; 19   : {

	movaps	xmm6, xmm3
	mov	edi, edx

; 50   : 	switch (mdly_modType)

	mov	DWORD PTR [rbx+108], ecx
	test	r8d, r8d
	je	SHORT $LN8@init
	sub	r8d, 1
	je	SHORT $LN7@init
	mov	DWORD PTR [rbx+112], 1056964608		; 3f000000H
	mov	DWORD PTR [rbx+92], 2
	cmp	r8d, 1
	jne	SHORT $LN14@init

; 51   : 	{
; 52   : 	case ModulationType::Flanger:
; 53   : 	{
; 54   : 		mdly_minDelaymSec = 1.0;
; 55   : 		mdly_maxDelaymSec = 7.0;
; 56   : 		mdly_wet = 0.5;
; 57   : 		mdly_dry = 0.5;
; 58   : 		mdly_feedback = 0.0;
; 59   : 		mdly_modLFO = OscillatorType::Triangular;	
; 60   : 		mdly_isUnipolar = true;	
; 61   : 		mdly_meanMod  = mdly_minDelaymSec;
; 62   : 		mdly_deltaMod = (mdly_maxDelaymSec - mdly_minDelaymSec);

	movss	xmm0, DWORD PTR __real@40c00000
	mov	DWORD PTR [rbx+100], 1065353216		; 3f800000H
	mov	DWORD PTR [rbx+104], 1088421888		; 40e00000H
	mov	DWORD PTR [rbx+116], 1056964608		; 3f000000H
	mov	BYTE PTR [rbx+88], r8b
	mov	DWORD PTR [rbx+120], 1065353216		; 3f800000H

; 63   : 		break;

	jmp	SHORT $LN4@init
$LN7@init:

; 64   : 	}
; 65   : 
; 66   : 	case ModulationType::Vibrato:
; 67   : 	{
; 68   : 		mdly_minDelaymSec = 0.0;
; 69   : 		mdly_maxDelaymSec = 7.0;
; 70   : 		mdly_wet = 1.0;
; 71   : 		mdly_dry = 0.0;
; 72   : 		mdly_feedback = 0.0;
; 73   : 		mdly_modLFO = OscillatorType::Sine;
; 74   : 		mdly_isUnipolar = true;
; 75   : 		mdly_meanMod = mdly_minDelaymSec;
; 76   : 		mdly_deltaMod = (mdly_maxDelaymSec - mdly_minDelaymSec);

	movss	xmm0, DWORD PTR __real@40e00000
	mov	DWORD PTR [rbx+100], ecx
	mov	DWORD PTR [rbx+104], 1088421888		; 40e00000H
	mov	QWORD PTR [rbx+112], 1065353216		; 3f800000H
	mov	DWORD PTR [rbx+92], 1
	mov	BYTE PTR [rbx+88], 1
	mov	DWORD PTR [rbx+120], ecx

; 77   : 		break;

	jmp	SHORT $LN4@init
$LN8@init:

; 78   : 	}
; 79   : 
; 80   : 	case ModulationType::Chorus:
; 81   : 	{
; 82   : 		mdly_minDelaymSec = 5.0;
; 83   : 		mdly_maxDelaymSec = 30.0;
; 84   : 		mdly_wet = 0.5;

	mov	DWORD PTR [rbx+112], 1056964608		; 3f000000H

; 85   : 		mdly_dry = 1.0;
; 86   : 		mdly_feedback = 0.0;
; 87   : 		mdly_modLFO = OscillatorType::Triangular;

	mov	DWORD PTR [rbx+92], 2
$LN14@init:

; 20   : 	// set modulation type
; 21   : 	setModType(modType);
; 22   : 
; 23   : 	// initialize delay line
; 24   : 	CombFilter::init(100.0, sampleRate);

	movss	xmm0, DWORD PTR __real@41480000
	mov	DWORD PTR [rbx+120], 1099694080		; 418c0000H
	mov	BYTE PTR [rbx+88], cl
	mov	DWORD PTR [rbx+116], 1065353216		; 3f800000H
	mov	DWORD PTR [rbx+104], 1106247680		; 41f00000H
	mov	DWORD PTR [rbx+100], 1084227584		; 40a00000H
$LN4@init:
	movss	DWORD PTR [rbx+124], xmm0
	mov	r8d, edi
	movss	xmm1, DWORD PTR __real@42c80000
	mov	rcx, rbx
	call	?init@Delay@@UEAAXMH@Z			; Delay::init
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\Delay.cpp

; 128  : 	dly_delayInmsec = delayInmsec;

	movss	xmm0, DWORD PTR [rbx+120]
	movd	xmm1, DWORD PTR [rbx+16]

; 129  : 
; 130  : 	if (dly_delayInmsec > dly_lineLengthInmsec)

	cvtdq2ps xmm1, xmm1
	movss	DWORD PTR [rbx+24], xmm0
	comiss	xmm0, xmm1
	jbe	SHORT $LN12@init

; 131  : 		dly_delayInmsec = dly_lineLengthInmsec;

	movss	DWORD PTR [rbx+24], xmm1
$LN12@init:

; 132  : 
; 133  : 	// Update parameters based on new delay length
; 134  : 	updateParameters();

	mov	rcx, rbx
	call	?updateParameters@Delay@@QEAAXXZ	; Delay::updateParameters
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp

; 30   : 	mdly_lfo->init(mdly_modLFO, mdly_modRate, sampleRate);

	mov	rcx, QWORD PTR [rbx+72]
	mov	r9d, edi
	mov	edx, DWORD PTR [rbx+92]
	movaps	xmm2, xmm6
	movss	DWORD PTR [rbx+80], xmm6
	mov	rax, QWORD PTR [rcx]
	call	QWORD PTR [rax]

; 31   : 	setModDepth(modDepth);
; 32   : 	mdly_lfo->setLFOunipolar(mdly_isUnipolar);

	mov	rcx, QWORD PTR [rbx+72]
	mov	al, BYTE PTR [rbx+88]
	movss	xmm0, DWORD PTR modDepth$[rsp]

; 33   : }

	movaps	xmm6, XMMWORD PTR [rsp+32]
	movss	DWORD PTR [rbx+84], xmm0
	mov	rbx, QWORD PTR [rsp+64]
	mov	BYTE PTR [rcx+64], al
	add	rsp, 48					; 00000030H
	pop	rdi
	ret	0
?init@ModDelay@@UEAAXHW4ModulationType@@MM@Z ENDP	; ModDelay::init
_TEXT	ENDS
; Function compile flags: /Ogspy
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp
;	COMDAT ?setModDepth@ModDelay@@QEAAXM@Z
_TEXT	SEGMENT
this$ = 8
modDepth$ = 16
?setModDepth@ModDelay@@QEAAXM@Z PROC			; ModDelay::setModDepth, COMDAT

; 43   : 	mdly_modDepth = modDepth;

	movss	DWORD PTR [rcx+84], xmm1

; 44   : 	//mdly_lfo->setLFOAmplitude(mdly_modDepth);
; 45   : }

	ret	0
?setModDepth@ModDelay@@QEAAXM@Z ENDP			; ModDelay::setModDepth
_TEXT	ENDS
; Function compile flags: /Ogspy
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp
;	COMDAT ?setModType@ModDelay@@QEAAXW4ModulationType@@@Z
_TEXT	SEGMENT
this$ = 8
modType$ = 16
?setModType@ModDelay@@QEAAXW4ModulationType@@@Z PROC	; ModDelay::setModType, COMDAT

; 49   : 	mdly_modType = modType;
; 50   : 	switch (mdly_modType)

	xor	r8d, r8d
	mov	DWORD PTR [rcx+96], edx
	mov	DWORD PTR [rcx+108], r8d
	test	edx, edx
	je	SHORT $LN6@setModType
	sub	edx, 1
	je	SHORT $LN5@setModType
	mov	DWORD PTR [rcx+112], 1056964608		; 3f000000H
	mov	DWORD PTR [rcx+92], 2
	cmp	edx, 1
	jne	SHORT $LN9@setModType

; 51   : 	{
; 52   : 	case ModulationType::Flanger:
; 53   : 	{
; 54   : 		mdly_minDelaymSec = 1.0;
; 55   : 		mdly_maxDelaymSec = 7.0;
; 56   : 		mdly_wet = 0.5;
; 57   : 		mdly_dry = 0.5;
; 58   : 		mdly_feedback = 0.0;
; 59   : 		mdly_modLFO = OscillatorType::Triangular;	
; 60   : 		mdly_isUnipolar = true;	
; 61   : 		mdly_meanMod  = mdly_minDelaymSec;
; 62   : 		mdly_deltaMod = (mdly_maxDelaymSec - mdly_minDelaymSec);

	movss	xmm0, DWORD PTR __real@40c00000
	mov	DWORD PTR [rcx+100], 1065353216		; 3f800000H
	mov	DWORD PTR [rcx+104], 1088421888		; 40e00000H
	mov	DWORD PTR [rcx+116], 1056964608		; 3f000000H
	mov	BYTE PTR [rcx+88], dl
	mov	DWORD PTR [rcx+120], 1065353216		; 3f800000H

; 63   : 		break;

	jmp	SHORT $LN2@setModType
$LN5@setModType:

; 64   : 	}
; 65   : 
; 66   : 	case ModulationType::Vibrato:
; 67   : 	{
; 68   : 		mdly_minDelaymSec = 0.0;
; 69   : 		mdly_maxDelaymSec = 7.0;
; 70   : 		mdly_wet = 1.0;
; 71   : 		mdly_dry = 0.0;
; 72   : 		mdly_feedback = 0.0;
; 73   : 		mdly_modLFO = OscillatorType::Sine;
; 74   : 		mdly_isUnipolar = true;
; 75   : 		mdly_meanMod = mdly_minDelaymSec;
; 76   : 		mdly_deltaMod = (mdly_maxDelaymSec - mdly_minDelaymSec);

	movss	xmm0, DWORD PTR __real@40e00000
	mov	DWORD PTR [rcx+100], r8d
	mov	DWORD PTR [rcx+104], 1088421888		; 40e00000H
	mov	QWORD PTR [rcx+112], 1065353216		; 3f800000H
	mov	DWORD PTR [rcx+92], 1
	mov	BYTE PTR [rcx+88], 1
	mov	DWORD PTR [rcx+120], r8d

; 77   : 		break;

	jmp	SHORT $LN2@setModType
$LN6@setModType:

; 78   : 	}
; 79   : 
; 80   : 	case ModulationType::Chorus:
; 81   : 	{
; 82   : 		mdly_minDelaymSec = 5.0;
; 83   : 		mdly_maxDelaymSec = 30.0;
; 84   : 		mdly_wet = 0.5;

	mov	DWORD PTR [rcx+112], 1056964608		; 3f000000H

; 85   : 		mdly_dry = 1.0;
; 86   : 		mdly_feedback = 0.0;
; 87   : 		mdly_modLFO = OscillatorType::Triangular;

	mov	DWORD PTR [rcx+92], 2
$LN9@setModType:

; 88   : 		mdly_isUnipolar = false;
; 89   : 		mdly_meanMod = (mdly_maxDelaymSec + mdly_minDelaymSec) / 2.0;
; 90   : 		mdly_deltaMod = (mdly_maxDelaymSec - mdly_minDelaymSec) / 2.0;
; 91   : 		break;
; 92   : 	}
; 93   : 
; 94   : 	default: // Chorus
; 95   : 	{
; 96   : 		mdly_minDelaymSec = 5.0;
; 97   : 		mdly_maxDelaymSec = 30.0;
; 98   : 		mdly_wet = 0.5;
; 99   : 		mdly_dry = 1.0;
; 100  : 		mdly_feedback = 0.0;
; 101  : 		mdly_modLFO = OscillatorType::Triangular;
; 102  : 		mdly_isUnipolar = false;
; 103  : 		mdly_meanMod = (mdly_maxDelaymSec + mdly_minDelaymSec) / 2.0;
; 104  : 		mdly_deltaMod = (mdly_maxDelaymSec - mdly_minDelaymSec) / 2.0;
; 105  : 		break;
; 106  : 	}
; 107  : 	}			
; 108  : }

	movss	xmm0, DWORD PTR __real@41480000
	mov	DWORD PTR [rcx+120], 1099694080		; 418c0000H
	mov	BYTE PTR [rcx+88], r8b
	mov	DWORD PTR [rcx+116], 1065353216		; 3f800000H
	mov	DWORD PTR [rcx+104], 1106247680		; 41f00000H
	mov	DWORD PTR [rcx+100], 1084227584		; 40a00000H
$LN2@setModType:
	movss	DWORD PTR [rcx+124], xmm0
	ret	0
?setModType@ModDelay@@QEAAXW4ModulationType@@@Z ENDP	; ModDelay::setModType
_TEXT	ENDS
; Function compile flags: /Ogspy
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp
;	COMDAT ?setSampleRate@ModDelay@@UEAAXH@Z
_TEXT	SEGMENT
this$ = 48
sampleRate$ = 56
?setSampleRate@ModDelay@@UEAAXH@Z PROC			; ModDelay::setSampleRate, COMDAT

; 111  : {

$LN4:
	mov	QWORD PTR [rsp+8], rbx
	push	rdi
	sub	rsp, 32					; 00000020H
	mov	rdi, rcx
	mov	ebx, edx

; 112  : 	mdly_lfo->setSampleRate(sampleRate);

	mov	rcx, QWORD PTR [rcx+72]
	mov	rax, QWORD PTR [rcx]
	call	QWORD PTR [rax+8]

; 113  : 	CombFilter::setSampleRate(sampleRate);

	mov	edx, ebx
	mov	rcx, rdi

; 114  : }

	mov	rbx, QWORD PTR [rsp+48]
	add	rsp, 32					; 00000020H
	pop	rdi

; 113  : 	CombFilter::setSampleRate(sampleRate);

	jmp	?setSampleRate@Delay@@UEAAXH@Z		; Delay::setSampleRate
?setSampleRate@ModDelay@@UEAAXH@Z ENDP			; ModDelay::setSampleRate
_TEXT	ENDS
; Function compile flags: /Ogspy
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\Delay.cpp
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\Delay.cpp
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp
;	COMDAT ?processAudio@ModDelay@@UEAAMM@Z
_TEXT	SEGMENT
this$ = 64
xn$ = 72
?processAudio@ModDelay@@UEAAMM@Z PROC			; ModDelay::processAudio, COMDAT

; 122  : {

$LN7:
	push	rbx
	sub	rsp, 48					; 00000030H
	mov	rbx, rcx
	movaps	XMMWORD PTR [rsp+32], xmm6

; 124  : 	float newDelayInmsec = mdly_meanMod +mdly_modDepth * mdly_deltaMod * mdly_lfo->processAudio();

	mov	rcx, QWORD PTR [rcx+72]
	movaps	xmm6, xmm1
	mov	rax, QWORD PTR [rcx]
	call	QWORD PTR [rax+16]
	movss	xmm2, DWORD PTR [rbx+124]
	mulss	xmm2, DWORD PTR [rbx+84]
	mulss	xmm2, xmm0
	movd	xmm0, DWORD PTR [rbx+16]
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\Delay.cpp

; 130  : 	if (dly_delayInmsec > dly_lineLengthInmsec)

	cvtdq2ps xmm0, xmm0
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp

; 124  : 	float newDelayInmsec = mdly_meanMod +mdly_modDepth * mdly_deltaMod * mdly_lfo->processAudio();

	addss	xmm2, DWORD PTR [rbx+120]
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\Delay.cpp

; 130  : 	if (dly_delayInmsec > dly_lineLengthInmsec)

	comiss	xmm2, xmm0
	movss	DWORD PTR [rbx+24], xmm2
	jbe	SHORT $LN4@processAud

; 131  : 		dly_delayInmsec = dly_lineLengthInmsec;

	movss	DWORD PTR [rbx+24], xmm0
$LN4@processAud:

; 132  : 
; 133  : 	// Update parameters based on new delay length
; 134  : 	updateParameters();

	mov	rcx, rbx
	call	?updateParameters@Delay@@QEAAXXZ	; Delay::updateParameters
; File C:\Users\RenatoVolpe\personal\fox-suite-blocks\src\ModDelay.cpp

; 130  : 	float yn = CombFilter::processAudio(xn);

	movaps	xmm1, xmm6
	mov	rcx, rbx
	call	?processAudio@CombFilter@@UEAAMM@Z	; CombFilter::processAudio

; 131  : 
; 132  : 	//// read from delay line
; 133  : 	//float buf = readFromDelayLine();
; 134  : 	//
; 135  : 	//// write do delay line
; 136  : 	//writeToDelayLine(mdly_feedback * buf + xn);
; 137  : 
; 138  : 	//// update delay indices
; 139  : 	//updateIndices();
; 140  : 
; 141  : 	return mdly_dry * xn + mdly_wet * yn;

	mulss	xmm6, DWORD PTR [rbx+116]
	mulss	xmm0, DWORD PTR [rbx+112]
	addss	xmm0, xmm6

; 142  : }

	movaps	xmm6, XMMWORD PTR [rsp+32]
	add	rsp, 48					; 00000030H
	pop	rbx
	ret	0
?processAudio@ModDelay@@UEAAMM@Z ENDP			; ModDelay::processAudio
_TEXT	ENDS
END
