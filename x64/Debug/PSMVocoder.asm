; Listing generated by Microsoft (R) Optimizing Compiler Version 19.29.30133.0 

include listing.inc

INCLUDELIB OLDNAMES

EXTRN	__imp_fftw_free:PROC
EXTRN	__imp_memset:PROC
EXTRN	??_U@YAPEAX_K@Z:PROC				; operator new[]
EXTRN	__imp_fftw_execute:PROC
EXTRN	__imp_fftw_destroy_plan:PROC
EXTRN	__imp_fftw_plan_dft_1d:PROC
EXTRN	__imp_fftw_malloc:PROC
?_Fake_alloc@std@@3U_Fake_allocator@1@B	ORG $+1		; std::_Fake_alloc
PUBLIC	?doOverlapAdd@PhaseVocoder@@QEAAXPEANH@Z	; PhaseVocoder::doOverlapAdd
PUBLIC	?doInverseFFT@PhaseVocoder@@QEAAXXZ		; PhaseVocoder::doInverseFFT
PUBLIC	?processAudioSample@PhaseVocoder@@QEAANNAEA_N@Z	; PhaseVocoder::processAudioSample
PUBLIC	?advanceAndCheckFFT@PhaseVocoder@@QEAA_NXZ	; PhaseVocoder::advanceAndCheckFFT
PUBLIC	?initialize@PhaseVocoder@@QEAAXIIW4windowType@@@Z ; PhaseVocoder::initialize
PUBLIC	?destroyFFTW@PhaseVocoder@@QEAAXXZ		; PhaseVocoder::destroyFFTW
PUBLIC	??_H@YAXPEAX_K1P6APEAX0@Z@Z			; `vector constructor iterator'
?kCTCorrFactorAntiLog@@3NB DQ 01H DUP (?)		; kCTCorrFactorAntiLog
?kCTCorrFactorZero@@3NB DQ 01H DUP (?)			; kCTCorrFactorZero
_BSS	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$?doInverseFFT@PhaseVocoder@@QEAAXXZ DD imagerel $LN4
	DD	imagerel $LN4+29
	DD	imagerel $unwind$?doInverseFFT@PhaseVocoder@@QEAAXXZ
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$?processAudioSample@PhaseVocoder@@QEAANNAEA_N@Z DD imagerel $LN8
	DD	imagerel $LN8+160
	DD	imagerel $unwind$?processAudioSample@PhaseVocoder@@QEAANNAEA_N@Z
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$?advanceAndCheckFFT@PhaseVocoder@@QEAA_NXZ DD imagerel $LN18
	DD	imagerel $LN18+176
	DD	imagerel $unwind$?advanceAndCheckFFT@PhaseVocoder@@QEAA_NXZ
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$?initialize@PhaseVocoder@@QEAAXIIW4windowType@@@Z DD imagerel $LN64
	DD	imagerel $LN64+1191
	DD	imagerel $unwind$?initialize@PhaseVocoder@@QEAAXIIW4windowType@@@Z
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$?destroyFFTW@PhaseVocoder@@QEAAXXZ DD imagerel $LN9
	DD	imagerel $LN9+89
	DD	imagerel $unwind$?destroyFFTW@PhaseVocoder@@QEAAXXZ
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$??_H@YAXPEAX_K1P6APEAX0@Z@Z DD imagerel $LN13
	DD	imagerel $LN13+47
	DD	imagerel $unwind$??_H@YAXPEAX_K1P6APEAX0@Z@Z
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$??__EkCTCorrFactorAntiLog@@YAXXZ DD imagerel ??__EkCTCorrFactorAntiLog@@YAXXZ
	DD	imagerel ??__EkCTCorrFactorAntiLog@@YAXXZ+47
	DD	imagerel $unwind$??__EkCTCorrFactorAntiLog@@YAXXZ
pdata	ENDS
;	COMDAT pdata
pdata	SEGMENT
$pdata$??__EkCTCorrFactorZero@@YAXXZ DD imagerel ??__EkCTCorrFactorZero@@YAXXZ
	DD	imagerel ??__EkCTCorrFactorZero@@YAXXZ+39
	DD	imagerel $unwind$??__EkCTCorrFactorZero@@YAXXZ
;	COMDAT xdata
xdata	SEGMENT
$unwind$??__EkCTCorrFactorZero@@YAXXZ DD 010401H
	DD	04204H
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$??__EkCTCorrFactorAntiLog@@YAXXZ DD 010401H
	DD	04204H
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$??_H@YAXPEAX_K1P6APEAX0@Z@Z DD 040a01H
	DD	06340aH
	DD	07006320aH
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$?destroyFFTW@PhaseVocoder@@QEAAXXZ DD 020601H
	DD	030023206H
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$?initialize@PhaseVocoder@@QEAAXIIW4windowType@@@Z DD 0e4301H
	DD	037843H
	DD	04681fH
	DD	0106418H
	DD	0f5418H
	DD	0e3418H
	DD	0f0149218H
	DD	07010d012H
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$?advanceAndCheckFFT@PhaseVocoder@@QEAA_NXZ DD 020601H
	DD	030023206H
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$?processAudioSample@PhaseVocoder@@QEAANNAEA_N@Z DD 081e01H
	DD	02781eH
	DD	036816H
	DD	0a340aH
	DD	07006720aH
xdata	ENDS
;	COMDAT xdata
xdata	SEGMENT
$unwind$?doInverseFFT@PhaseVocoder@@QEAAXXZ DD 020601H
	DD	030023206H
?kCTCorrFactorZero$initializer$@@3P6AXXZEA DQ FLAT:??__EkCTCorrFactorZero@@YAXXZ ; kCTCorrFactorZero$initializer$
?kCTCorrFactorUnity$initializer$@@3P6AXXZEA DQ FLAT:??__EkCTCorrFactorUnity@@YAXXZ ; kCTCorrFactorUnity$initializer$
?kCTCorrFactorAntiUnity$initializer$@@3P6AXXZEA DQ FLAT:??__EkCTCorrFactorAntiUnity@@YAXXZ ; kCTCorrFactorAntiUnity$initializer$
?kCTCorrFactorAntiLog$initializer$@@3P6AXXZEA DQ FLAT:??__EkCTCorrFactorAntiLog@@YAXXZ ; kCTCorrFactorAntiLog$initializer$
?kCTCorrFactorAntiLogScale$initializer$@@3P6AXXZEA DQ FLAT:??__EkCTCorrFactorAntiLogScale@@YAXXZ ; kCTCorrFactorAntiLogScale$initializer$
; Function compile flags: /Ogspy
; File E:\prova\Shimmer\include\guiconstants.h
;	COMDAT ??__EkCTCorrFactorZero@@YAXXZ
text$di	SEGMENT
??__EkCTCorrFactorZero@@YAXXZ PROC			; `dynamic initializer for 'kCTCorrFactorZero'', COMDAT

; 142  : const double kCTCorrFactorZero = pow(10.0, (-1.0/kCTCoefficient));

	sub	rsp, 40					; 00000028H
	movsd	xmm1, QWORD PTR __real@c003333333333333
	movsd	xmm0, QWORD PTR __real@4024000000000000
	call	QWORD PTR __imp_pow
	movsd	QWORD PTR ?kCTCorrFactorZero@@3NB, xmm0
	add	rsp, 40					; 00000028H
	ret	0
??__EkCTCorrFactorZero@@YAXXZ ENDP			; `dynamic initializer for 'kCTCorrFactorZero''
text$di	ENDS
; Function compile flags: /Ogspy
; File E:\prova\Shimmer\include\guiconstants.h
;	COMDAT ??__EkCTCorrFactorUnity@@YAXXZ
text$di	SEGMENT
??__EkCTCorrFactorUnity@@YAXXZ PROC			; `dynamic initializer for 'kCTCorrFactorUnity'', COMDAT

; 156  : const double kCTCorrFactorUnity = 1.0 / (1.0 + kCTCoefficient*log10(1.0 + kCTCorrFactorZero));

	movsd	xmm0, QWORD PTR ?kCTCorrFactorZero@@3NB
	addsd	xmm0, QWORD PTR __real@3ff0000000000000
	rex_jmp	QWORD PTR __imp_log10
??__EkCTCorrFactorUnity@@YAXXZ ENDP			; `dynamic initializer for 'kCTCorrFactorUnity''
text$di	ENDS
; Function compile flags: /Ogspy
; File E:\prova\Shimmer\include\guiconstants.h
;	COMDAT ??__EkCTCorrFactorAntiUnity@@YAXXZ
text$di	SEGMENT
??__EkCTCorrFactorAntiUnity@@YAXXZ PROC			; `dynamic initializer for 'kCTCorrFactorAntiUnity'', COMDAT

; 163  : const double kCTCorrFactorAntiUnity = 1.0 / (1.0 + (-pow(10.0, (-1.0/kCTCoefficient))));

	movsd	xmm1, QWORD PTR __real@c003333333333333
	movsd	xmm0, QWORD PTR __real@4024000000000000
	rex_jmp	QWORD PTR __imp_pow
??__EkCTCorrFactorAntiUnity@@YAXXZ ENDP			; `dynamic initializer for 'kCTCorrFactorAntiUnity''
text$di	ENDS
; Function compile flags: /Ogspy
; File E:\prova\Shimmer\include\guiconstants.h
;	COMDAT ??__EkCTCorrFactorAntiLog@@YAXXZ
text$di	SEGMENT
??__EkCTCorrFactorAntiLog@@YAXXZ PROC			; `dynamic initializer for 'kCTCorrFactorAntiLog'', COMDAT

; 170  : const double kCTCorrFactorAntiLog = kCTCoefficient*log10(1.0 + kCTCorrFactorZero);

	sub	rsp, 40					; 00000028H
	movsd	xmm0, QWORD PTR ?kCTCorrFactorZero@@3NB
	addsd	xmm0, QWORD PTR __real@3ff0000000000000
	call	QWORD PTR __imp_log10
	mulsd	xmm0, QWORD PTR __real@3fdaaaaaaaaaaaab
	movsd	QWORD PTR ?kCTCorrFactorAntiLog@@3NB, xmm0
	add	rsp, 40					; 00000028H
	ret	0
??__EkCTCorrFactorAntiLog@@YAXXZ ENDP			; `dynamic initializer for 'kCTCorrFactorAntiLog''
text$di	ENDS
; Function compile flags: /Ogspy
; File E:\prova\Shimmer\include\guiconstants.h
;	COMDAT ??__EkCTCorrFactorAntiLogScale@@YAXXZ
text$di	SEGMENT
??__EkCTCorrFactorAntiLogScale@@YAXXZ PROC		; `dynamic initializer for 'kCTCorrFactorAntiLogScale'', COMDAT

; 177  : const double kCTCorrFactorAntiLogScale = 1.0 / (-kCTCoefficient*log10(kCTCorrFactorZero) + kCTCorrFactorAntiLog);

	movsd	xmm0, QWORD PTR ?kCTCorrFactorZero@@3NB
	rex_jmp	QWORD PTR __imp_log10
??__EkCTCorrFactorAntiLogScale@@YAXXZ ENDP		; `dynamic initializer for 'kCTCorrFactorAntiLogScale''
text$di	ENDS
; Function compile flags: /Ogspy
;	COMDAT ??_H@YAXPEAX_K1P6APEAX0@Z@Z
_TEXT	SEGMENT
__t$ = 48
__s$dead$ = 56
__n$dead$ = 64
__f$dead$ = 72
??_H@YAXPEAX_K1P6APEAX0@Z@Z PROC			; `vector constructor iterator', COMDAT
$LN13:
	mov	QWORD PTR [rsp+8], rbx
	push	rdi
	sub	rsp, 32					; 00000020H
	mov	rbx, rcx
	mov	edi, 4096				; 00001000H
$LL2@vector:
	mov	rcx, rbx
	call	??0BinData@@QEAA@XZ			; BinData::BinData
	add	rbx, 48					; 00000030H
	sub	rdi, 1
	jne	SHORT $LL2@vector
	mov	rbx, QWORD PTR [rsp+48]
	add	rsp, 32					; 00000020H
	pop	rdi
	ret	0
??_H@YAXPEAX_K1P6APEAX0@Z@Z ENDP			; `vector constructor iterator'
_TEXT	ENDS
; Function compile flags: /Ogspy
; File E:\prova\Shimmer\include\src\PSMVocoder.cpp
;	COMDAT ?destroyFFTW@PhaseVocoder@@QEAAXXZ
_TEXT	SEGMENT
this$ = 48
?destroyFFTW@PhaseVocoder@@QEAAXXZ PROC			; PhaseVocoder::destroyFFTW, COMDAT

; 186  : {

$LN9:
	push	rbx
	sub	rsp, 32					; 00000020H
	mov	rbx, rcx

; 187  : 	if (plan_forward)

	mov	rcx, QWORD PTR [rcx+24]
	test	rcx, rcx
	je	SHORT $LN2@destroyFFT

; 188  : 		fftw_destroy_plan(plan_forward);

	call	QWORD PTR __imp_fftw_destroy_plan
$LN2@destroyFFT:

; 189  : 	if (plan_backward)

	mov	rcx, QWORD PTR [rbx+32]
	test	rcx, rcx
	je	SHORT $LN3@destroyFFT

; 190  : 		fftw_destroy_plan(plan_backward);

	call	QWORD PTR __imp_fftw_destroy_plan
$LN3@destroyFFT:

; 191  : 
; 192  : 	if (fft_input)

	mov	rcx, QWORD PTR [rbx]
	test	rcx, rcx
	je	SHORT $LN4@destroyFFT

; 193  : 		fftw_free(fft_input);

	call	QWORD PTR __imp_fftw_free
$LN4@destroyFFT:

; 194  : 	if (fft_result)

	mov	rcx, QWORD PTR [rbx+8]
	test	rcx, rcx
	je	SHORT $LN5@destroyFFT

; 195  : 		fftw_free(fft_result);

	call	QWORD PTR __imp_fftw_free
$LN5@destroyFFT:

; 196  : 	if (ifft_result)

	mov	rcx, QWORD PTR [rbx+16]
	test	rcx, rcx
	je	SHORT $LN6@destroyFFT

; 197  : 		fftw_free(ifft_result);

	call	QWORD PTR __imp_fftw_free
$LN6@destroyFFT:

; 198  : }

	add	rsp, 32					; 00000020H
	pop	rbx
	ret	0
?destroyFFTW@PhaseVocoder@@QEAAXXZ ENDP			; PhaseVocoder::destroyFFTW
_TEXT	ENDS
; Function compile flags: /Ogspy
; File E:\prova\Shimmer\include\src\PSMVocoder.cpp
;	COMDAT ?initialize@PhaseVocoder@@QEAAXIIW4windowType@@@Z
_TEXT	SEGMENT
this$ = 112
_frameLength$dead$ = 120
_hopSize$dead$ = 128
_window$dead$ = 136
?initialize@PhaseVocoder@@QEAAXIIW4windowType@@@Z PROC	; PhaseVocoder::initialize, COMDAT

; 212  : {

$LN64:
	mov	rax, rsp
	mov	QWORD PTR [rax+8], rbx
	mov	QWORD PTR [rax+16], rbp
	mov	QWORD PTR [rax+24], rsi
	push	rdi
	push	r13
	push	r15
	sub	rsp, 80					; 00000050H
	mov	rbx, rcx
	movaps	XMMWORD PTR [rax-40], xmm6

; 213  : 	frameLength = _frameLength;
; 214  : 	wrapMask = frameLength - 1;

	mov	DWORD PTR [rcx+80], 4095		; 00000fffH

; 215  : 	hopSize = _hopSize;
; 216  : 	window = _window;
; 217  : 
; 218  : 	// --- this is the overlap as a fraction i.e. 0.75 = 75%
; 219  : 	overlap = hopSize > 0.0 ? 1.0 - (double)hopSize / (double)frameLength : 0.0;
; 220  : 
; 221  : 	// --- gain correction for window + hop size
; 222  : 	windowHopCorrection = 0.0;
; 223  : 
; 224  : 	// --- SETUP BUFFERS ---- //
; 225  : 	//     NOTE: input and output buffers are circular, others are linear
; 226  : 	//
; 227  : 	// --- input buffer, for processing the x(n) timeline
; 228  : 	if (inputBuffer)

	mov	r15d, 8
	mov	DWORD PTR [rcx+112], 1024		; 00000400H
	mov	DWORD PTR [rcx+100], 2
	and	QWORD PTR [rbx+88], 0
	movaps	XMMWORD PTR [rax-56], xmm7
	mov	eax, 4096				; 00001000H
	mov	DWORD PTR [rcx+104], eax
	mov	rcx, 4604930618986332160		; 3fe8000000000000H
	mov	QWORD PTR [rbx+120], rcx
	mov	rcx, QWORD PTR [rbx+48]
	test	rcx, rcx
	je	SHORT $LN20@initialize

; 229  : 		delete inputBuffer;

	mov	edx, r15d
	call	??3@YAXPEAX_K@Z				; operator delete
	mov	eax, DWORD PTR [rbx+104]
$LN20@initialize:

; 230  : 
; 231  : 	inputBuffer = new double[frameLength];

	mov	ecx, eax
	mov	r13, -1
	mov	rax, r15
	mul	rcx
	cmovo	rax, r13
	mov	rcx, rax
	call	??_U@YAPEAX_K@Z				; operator new[]

; 232  : 	memset(&inputBuffer[0], 0, frameLength * sizeof(double));

	mov	r8d, DWORD PTR [rbx+104]
	xor	edx, edx
	shl	r8, 3
	mov	rcx, rax
	mov	QWORD PTR [rbx+48], rax
	call	QWORD PTR __imp_memset

; 233  : 
; 234  : 	// --- output buffer, for processing the y(n) timeline and accumulating frames
; 235  : 	if (outputBuffer)

	mov	rcx, QWORD PTR [rbx+56]
	test	rcx, rcx
	je	SHORT $LN21@initialize

; 236  : 		delete outputBuffer;

	mov	rdx, r15
	call	??3@YAXPEAX_K@Z				; operator delete
$LN21@initialize:

; 237  : 
; 238  : 	// --- the output buffer is declared as 2x the normal frame size
; 239  : 	//     to accomodate time-stretching/pitch shifting; you can increase the size
; 240  : 	//     here; if so make sure to calculate the wrapMaskOut properly and everything
; 241  : 	//     will work normally you can even dynamically expand and contract the buffer
; 242  : 	//     (not sure why you would do this - and it will surely affect CPU performance)
; 243  : 	//     NOTE: the length of the buffer is only to accomodate accumulations
; 244  : 	//           it does not stretch time or change causality on its own
; 245  : 	outputBuffer = new double[frameLength * 4];

	mov	ecx, DWORD PTR [rbx+104]
	mov	rax, r15
	shl	ecx, 2
	mul	rcx
	cmovo	rax, r13
	mov	rcx, rax
	call	??_U@YAPEAX_K@Z				; operator new[]

; 246  : 	memset(&outputBuffer[0], 0, (frameLength * 4.0) * sizeof(double));

	mov	ecx, DWORD PTR [rbx+104]
	xorps	xmm0, xmm0
	movsd	xmm1, QWORD PTR __real@43e0000000000000
	mov	QWORD PTR [rbx+56], rax
	cvtsi2sd xmm0, rcx
	xor	ecx, ecx
	mulsd	xmm0, QWORD PTR __real@4040000000000000
	comisd	xmm0, xmm1
	jb	SHORT $LN62@initialize
	subsd	xmm0, xmm1
	comisd	xmm0, xmm1
	jae	SHORT $LN62@initialize
	mov	rdx, -9223372036854775808		; 8000000000000000H
	mov	rcx, rdx
$LN62@initialize:
	cvttsd2si r8, xmm0
	xor	edx, edx
	add	r8, rcx
	mov	rcx, rax
	call	QWORD PTR __imp_memset

; 247  : 	wrapMaskOut = (frameLength * 4.0) - 1;

	mov	edx, DWORD PTR [rbx+104]
	xorps	xmm0, xmm0
	movsd	xmm7, QWORD PTR __real@3ff0000000000000

; 248  : 
; 249  : 	// --- fixed window buffer
; 250  : 	if (windowBuffer)

	mov	rcx, QWORD PTR [rbx+40]
	cvtsi2sd xmm0, rdx
	mulsd	xmm0, QWORD PTR __real@4010000000000000
	subsd	xmm0, xmm7
	cvttsd2si rax, xmm0
	mov	DWORD PTR [rbx+84], eax
	test	rcx, rcx
	je	SHORT $LN22@initialize

; 251  : 		delete windowBuffer;

	mov	rdx, r15
	call	??3@YAXPEAX_K@Z				; operator delete
	mov	edx, DWORD PTR [rbx+104]
$LN22@initialize:

; 252  : 
; 253  : 	windowBuffer = new double[frameLength];

	mov	ecx, edx
	mov	rax, r15
	mul	rcx
	cmovo	rax, r13
	mov	rcx, rax
	call	??_U@YAPEAX_K@Z				; operator new[]

; 254  : 	memset(&windowBuffer[0], 0, frameLength * sizeof(double));

	mov	r8d, DWORD PTR [rbx+104]
	xor	edx, edx
	shl	r8, 3
	mov	rcx, rax
	mov	QWORD PTR [rbx+40], rax
	call	QWORD PTR __imp_memset

; 255  : 
; 256  : 	// --- this is from Reiss & McPherson's code
; 257  : 	//     https://code.soundsoftware.ac.uk/projects/audio_effects_textbook_code/repository/entry/effects/pvoc_passthrough/Source/PluginProcessor.cpp
; 258  : 	// NOTE:	"Window functions are typically defined to be symmetrical. This will cause a
; 259  : 	//			problem in the overlap-add process: the windows instead need to be periodic
; 260  : 	//			when arranged end-to-end. As a result we calculate the window of one sample
; 261  : 	//			larger than usual, and drop the last sample. (This works as long as N is even.)
; 262  : 	//			See Julius Smith, "Spectral Audio Signal Processing" for details.
; 263  : 	// --- WP: this is why denominators are (frameLength) rather than (frameLength - 1)
; 264  : 	if (window == windowType::kRectWindow)

	mov	ecx, DWORD PTR [rbx+100]
	mov	eax, DWORD PTR [rbx+104]
	cmp	ecx, 1
	jne	SHORT $LN23@initialize

; 265  : 	{
; 266  : 		for (int n = 0; n < frameLength - 1; n++)

	xor	r8d, r8d
	cmp	eax, ecx
	je	$LN18@initialize
	mov	r9, QWORD PTR [rbx+40]
	xor	edx, edx
	mov	rcx, 4607182418800017408		; 3ff0000000000000H
$LL4@initialize:

; 267  : 		{
; 268  : 			windowBuffer[n] = 1.0;

	mov	QWORD PTR [rdx+r9], rcx
	inc	r8d

; 269  : 			windowHopCorrection += windowBuffer[n];

	mov	r9, QWORD PTR [rbx+40]
	movsd	xmm0, QWORD PTR [rdx+r9]
	add	rdx, r15
	addsd	xmm0, QWORD PTR [rbx+88]
	movsd	QWORD PTR [rbx+88], xmm0
	mov	eax, DWORD PTR [rbx+104]
	dec	eax
	cmp	r8d, eax
	jb	SHORT $LL4@initialize

; 270  : 		}
; 271  : 	}

	jmp	$LN18@initialize
$LN23@initialize:

; 272  : 	else if (window == windowType::kHammingWindow)

	cmp	ecx, 4
	jne	SHORT $LN25@initialize

; 273  : 	{
; 274  : 		for (int n = 0; n < frameLength - 1; n++)

	xor	esi, esi
	cmp	eax, 1
	je	$LN18@initialize
	xor	edi, edi
$LL7@initialize:

; 275  : 		{
; 276  : 			windowBuffer[n] = 0.54 - 0.46 * cos((n * 2.0 * kPi) / (frameLength));

	mov	eax, DWORD PTR [rbx+104]
	xorps	xmm1, xmm1
	movd	xmm0, esi
	cvtdq2pd xmm0, xmm0
	cvtsi2sd xmm1, rax
	mulsd	xmm0, QWORD PTR __real@401921fb54442d18
	divsd	xmm0, xmm1
	call	QWORD PTR __imp_cos
	mov	rax, QWORD PTR [rbx+40]
	inc	esi
	mulsd	xmm0, QWORD PTR __real@3fdd70a3d70a3d71
	movsd	xmm1, QWORD PTR __real@3fe147ae147ae148
	subsd	xmm1, xmm0
	movsd	QWORD PTR [rdi+rax], xmm1

; 277  : 			windowHopCorrection += windowBuffer[n];

	mov	rax, QWORD PTR [rbx+40]
	movsd	xmm0, QWORD PTR [rdi+rax]
	add	rdi, r15
	addsd	xmm0, QWORD PTR [rbx+88]
	movsd	QWORD PTR [rbx+88], xmm0
	mov	eax, DWORD PTR [rbx+104]
	dec	eax
	cmp	esi, eax
	jb	SHORT $LL7@initialize

; 278  : 		}
; 279  : 	}

	jmp	$LN18@initialize
$LN25@initialize:

; 280  : 	else if (window == windowType::kHannWindow)

	cmp	ecx, 2
	jne	SHORT $LN27@initialize

; 281  : 	{
; 282  : 		for (int n = 0; n < frameLength; n++)

	xor	edi, edi
	test	eax, eax
	je	$LN18@initialize
	xor	esi, esi
$LL10@initialize:

; 283  : 		{
; 284  : 			windowBuffer[n] = 0.5 * (1 - cos((n * 2.0 * kPi) / (frameLength)));

	mov	eax, DWORD PTR [rbx+104]
	xorps	xmm1, xmm1
	movd	xmm0, edi
	cvtdq2pd xmm0, xmm0
	cvtsi2sd xmm1, rax
	mulsd	xmm0, QWORD PTR __real@401921fb54442d18
	divsd	xmm0, xmm1
	call	QWORD PTR __imp_cos
	mov	rax, QWORD PTR [rbx+40]
	movaps	xmm1, xmm7
	subsd	xmm1, xmm0
	inc	edi
	mulsd	xmm1, QWORD PTR __real@3fe0000000000000
	movsd	QWORD PTR [rsi+rax], xmm1

; 285  : 			windowHopCorrection += windowBuffer[n];

	mov	rax, QWORD PTR [rbx+40]
	movsd	xmm0, QWORD PTR [rsi+rax]
	add	rsi, r15
	addsd	xmm0, QWORD PTR [rbx+88]
	movsd	QWORD PTR [rbx+88], xmm0
	cmp	edi, DWORD PTR [rbx+104]
	jb	SHORT $LL10@initialize

; 286  : 		}
; 287  : 	}

	jmp	$LN18@initialize
$LN27@initialize:

; 288  : 	else if (window == windowType::kBlackmanHarrisWindow)

	cmp	ecx, 3
	jne	$LN29@initialize

; 289  : 	{
; 290  : 		for (int n = 0; n < frameLength; n++)

	xor	edi, edi
	test	eax, eax
	je	$LN18@initialize
	xor	esi, esi
	xor	ebp, ebp
$LL13@initialize:

; 291  : 		{
; 292  : 			windowBuffer[n] = (0.42323 - (0.49755 * cos((n * 2.0 * kPi) / (frameLength))) + 0.07922 * cos((2 * n * 2.0 * kPi) / (frameLength)));

	mov	eax, DWORD PTR [rbx+104]
	xorps	xmm0, xmm0
	cvtsi2sd xmm0, edi
	xorps	xmm1, xmm1
	cvtsi2sd xmm1, rax
	mulsd	xmm0, QWORD PTR __real@401921fb54442d18
	divsd	xmm0, xmm1
	call	QWORD PTR __imp_cos
	mov	eax, DWORD PTR [rbx+104]
	xorps	xmm1, xmm1
	mulsd	xmm0, QWORD PTR __real@3fdfd7dbf487fcb9
	movsd	xmm6, QWORD PTR __real@3fdb1633482be8bc
	cvtsi2sd xmm1, rax
	subsd	xmm6, xmm0
	xorps	xmm0, xmm0
	cvtsi2sd xmm0, ebp
	mulsd	xmm0, QWORD PTR __real@401921fb54442d18
	divsd	xmm0, xmm1
	call	QWORD PTR __imp_cos
	mov	rax, QWORD PTR [rbx+40]
	inc	edi
	mulsd	xmm0, QWORD PTR __real@3fb447c30d306a2b
	add	ebp, 2
	addsd	xmm0, xmm6
	movsd	QWORD PTR [rsi+rax], xmm0

; 293  : 			windowHopCorrection += windowBuffer[n];

	mov	rax, QWORD PTR [rbx+40]
	movsd	xmm0, QWORD PTR [rsi+rax]
	add	rsi, r15
	addsd	xmm0, QWORD PTR [rbx+88]
	movsd	QWORD PTR [rbx+88], xmm0
	cmp	edi, DWORD PTR [rbx+104]
	jb	$LL13@initialize

; 294  : 		}
; 295  : 	}

	jmp	SHORT $LN18@initialize
$LN29@initialize:

; 296  : 	else if (window == windowType::kNoWindow)

	xor	edx, edx
	test	ecx, ecx
	jne	SHORT $LN31@initialize

; 297  : 	{
; 298  : 		for (int n = 0; n < frameLength; n++)

	test	eax, eax
	je	SHORT $LN18@initialize
	mov	r8, QWORD PTR [rbx+40]
	xor	eax, eax
	mov	rcx, 4607182418800017408		; 3ff0000000000000H
$LL16@initialize:

; 299  : 		{
; 300  : 			windowBuffer[n] = 1.0;

	mov	QWORD PTR [rax+r8], rcx
	inc	edx

; 301  : 			windowHopCorrection += windowBuffer[n];

	mov	r8, QWORD PTR [rbx+40]
	movsd	xmm0, QWORD PTR [rax+r8]
	add	rax, r15
	addsd	xmm0, QWORD PTR [rbx+88]
	movsd	QWORD PTR [rbx+88], xmm0
	cmp	edx, DWORD PTR [rbx+104]
	jb	SHORT $LL16@initialize

; 302  : 		}
; 303  : 	}

	jmp	SHORT $LN18@initialize
$LN31@initialize:

; 304  : 	else // --- default to kNoWindow
; 305  : 	{
; 306  : 		for (int n = 0; n < frameLength; n++)

	test	eax, eax
	je	SHORT $LN18@initialize
	mov	r8, QWORD PTR [rbx+40]
	xor	eax, eax
	mov	rcx, 4607182418800017408		; 3ff0000000000000H
$LL19@initialize:

; 307  : 		{
; 308  : 			windowBuffer[n] = 1.0;

	mov	QWORD PTR [rax+r8], rcx
	inc	edx

; 309  : 			windowHopCorrection += windowBuffer[n];

	mov	r8, QWORD PTR [rbx+40]
	movsd	xmm0, QWORD PTR [rax+r8]
	add	rax, r15
	addsd	xmm0, QWORD PTR [rbx+88]
	movsd	QWORD PTR [rbx+88], xmm0
	cmp	edx, DWORD PTR [rbx+104]
	jb	SHORT $LL19@initialize
$LN18@initialize:

; 310  : 		}
; 311  : 	}
; 312  : 
; 313  : 	// --- calculate gain correction factor
; 314  : 	if (window != windowType::kNoWindow)

	cmp	DWORD PTR [rbx+100], 0
	je	SHORT $LN33@initialize

; 315  : 		windowHopCorrection = (1.0 - overlap) / windowHopCorrection;

	subsd	xmm7, QWORD PTR [rbx+120]
$LN33@initialize:

; 316  : 	else
; 317  : 		windowHopCorrection = 1.0 / windowHopCorrection;
; 318  : 
; 319  : 	// --- set
; 320  : 	inputWriteIndex = 0;

	divsd	xmm7, QWORD PTR [rbx+88]

; 321  : 	inputReadIndex = 0;
; 322  : 
; 323  : 	outputWriteIndex = 0;
; 324  : 	outputReadIndex = 0;
; 325  : 
; 326  : 	fftCounter = 0;
; 327  : 
; 328  : 	// --- reset flags
; 329  : 	needInverseFFT = false;
; 330  : 	needOverlapAdd = false;
; 331  : 
; 332  : #ifdef HAVE_FFTW
; 333  : 	destroyFFTW();

	mov	rcx, rbx
	movsd	QWORD PTR [rbx+88], xmm7
	and	DWORD PTR [rbx+64], 0
	and	DWORD PTR [rbx+72], 0
	and	DWORD PTR [rbx+68], 0
	and	DWORD PTR [rbx+76], 0
	and	DWORD PTR [rbx+108], 0
	and	WORD PTR [rbx+96], 0
	call	?destroyFFTW@PhaseVocoder@@QEAAXXZ	; PhaseVocoder::destroyFFTW

; 334  : 	fft_input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * frameLength);

	mov	ecx, DWORD PTR [rbx+104]
	shl	rcx, 4
	call	QWORD PTR __imp_fftw_malloc

; 335  : 	fft_result = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * frameLength);

	mov	ecx, DWORD PTR [rbx+104]
	shl	rcx, 4
	mov	QWORD PTR [rbx], rax
	call	QWORD PTR __imp_fftw_malloc

; 336  : 	ifft_result = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * frameLength);

	mov	ecx, DWORD PTR [rbx+104]
	shl	rcx, 4
	mov	QWORD PTR [rbx+8], rax
	call	QWORD PTR __imp_fftw_malloc

; 337  : 
; 338  : 	plan_forward = fftw_plan_dft_1d(frameLength, fft_input, fft_result, FFTW_FORWARD, FFTW_ESTIMATE);

	mov	r8, QWORD PTR [rbx+8]
	mov	edi, 64					; 00000040H
	mov	rdx, QWORD PTR [rbx]
	mov	r9d, r13d
	mov	ecx, DWORD PTR [rbx+104]
	mov	QWORD PTR [rbx+16], rax
	mov	DWORD PTR [rsp+32], edi
	call	QWORD PTR __imp_fftw_plan_dft_1d

; 339  : 	plan_backward = fftw_plan_dft_1d(frameLength, fft_result, ifft_result, FFTW_BACKWARD, FFTW_ESTIMATE);

	mov	r8, QWORD PTR [rbx+16]
	lea	r9d, QWORD PTR [rdi-63]
	mov	rdx, QWORD PTR [rbx+8]
	mov	ecx, DWORD PTR [rbx+104]
	mov	QWORD PTR [rbx+24], rax
	mov	DWORD PTR [rsp+32], edi
	call	QWORD PTR __imp_fftw_plan_dft_1d

; 340  : #endif
; 341  : }

	movaps	xmm6, XMMWORD PTR [rsp+64]
	lea	r11, QWORD PTR [rsp+80]
	mov	rbp, QWORD PTR [r11+40]
	mov	rsi, QWORD PTR [r11+48]
	movaps	xmm7, XMMWORD PTR [rsp+48]
	mov	QWORD PTR [rbx+32], rax
	mov	rbx, QWORD PTR [r11+32]
	mov	rsp, r11
	pop	r15
	pop	r13
	pop	rdi
	ret	0
?initialize@PhaseVocoder@@QEAAXIIW4windowType@@@Z ENDP	; PhaseVocoder::initialize
_TEXT	ENDS
; Function compile flags: /Ogspy
; File E:\prova\Shimmer\include\src\PSMVocoder.cpp
;	COMDAT ?advanceAndCheckFFT@PhaseVocoder@@QEAA_NXZ
_TEXT	SEGMENT
this$ = 48
?advanceAndCheckFFT@PhaseVocoder@@QEAA_NXZ PROC		; PhaseVocoder::advanceAndCheckFFT, COMDAT

; 383  : {

$LN18:
	push	rbx
	sub	rsp, 32					; 00000020H

; 384  : 	// --- inc counter and check count
; 385  : 	fftCounter++;

	inc	DWORD PTR [rcx+108]
	mov	rbx, rcx

; 386  : 
; 387  : 	if (fftCounter != frameLength)

	mov	eax, DWORD PTR [rcx+104]
	cmp	DWORD PTR [rcx+108], eax
	je	SHORT $LN5@advanceAnd

; 388  : 		return false;

	xor	al, al
	jmp	$LN1@advanceAnd
$LN5@advanceAnd:

; 389  : 
; 390  : 	// --- we have a FFT ready
; 391  : 	// --- load up the input to the FFT
; 392  : 	for (int i = 0; i < frameLength; i++)

	xor	r9d, r9d
	test	eax, eax
	je	SHORT $LN16@advanceAnd
	xor	r8d, r8d
	xor	r10d, r10d
$LL4@advanceAnd:

; 393  : 	{
; 394  : 		fft_input[i][0] = inputBuffer[inputReadIndex++] * windowBuffer[i];

	mov	edx, DWORD PTR [rbx+72]
	inc	r9d
	mov	rax, QWORD PTR [rbx+48]
	mov	rcx, QWORD PTR [rbx+40]
	movsd	xmm0, QWORD PTR [rax+rdx*8]
	mulsd	xmm0, QWORD PTR [rcx+r10]
	lea	r10, QWORD PTR [r10+8]
	mov	rax, QWORD PTR [rbx]
	movsd	QWORD PTR [r8+rax], xmm0
	lea	r8, QWORD PTR [r8+16]
	inc	DWORD PTR [rbx+72]

; 395  : 		fft_input[i][1] = 0.0; // use this if your data is complex valued

	mov	rax, QWORD PTR [rbx]
	and	QWORD PTR [r8+rax-8], 0

; 396  : 
; 397  : 		// --- wrap if index > bufferlength - 1
; 398  : 		inputReadIndex &= wrapMask;

	mov	eax, DWORD PTR [rbx+80]
	and	DWORD PTR [rbx+72], eax
	cmp	r9d, DWORD PTR [rbx+104]
	jb	SHORT $LL4@advanceAnd
$LN16@advanceAnd:

; 399  : 	}
; 400  : 
; 401  : 	// --- do the FFT
; 402  : 	fftw_execute(plan_forward);

	mov	rcx, QWORD PTR [rbx+24]
	call	QWORD PTR __imp_fftw_execute

; 403  : 
; 404  : 	// --- in case user does not take IFFT, just to prevent zero output
; 405  : 	needInverseFFT = true;
; 406  : 	needOverlapAdd = true;
; 407  : 
; 408  : 	// --- fft counter: small hop = more FFTs = less counting before fft
; 409  : 	//
; 410  : 	// --- overlap-add-only algorithms do not involve hop-size in FFT count
; 411  : 	if (overlapAddOnly)

	mov	cl, BYTE PTR [rbx+128]
	mov	WORD PTR [rbx+96], 257			; 00000101H
	test	cl, cl
	je	SHORT $LN6@advanceAnd

; 412  : 		fftCounter = 0;

	xor	eax, eax
	jmp	SHORT $LN7@advanceAnd
$LN6@advanceAnd:

; 413  : 	else // normal counter advance
; 414  : 		fftCounter = frameLength - hopSize;

	mov	eax, DWORD PTR [rbx+104]
	sub	eax, DWORD PTR [rbx+112]
$LN7@advanceAnd:

; 415  : 
; 416  : 	// --- setup the read index for next time through the loop
; 417  : 	if (!overlapAddOnly)

	mov	DWORD PTR [rbx+108], eax
	test	cl, cl
	jne	SHORT $LN14@advanceAnd

; 418  : 		inputReadIndex += hopSize;

	mov	eax, DWORD PTR [rbx+112]
	add	eax, DWORD PTR [rbx+72]
	jmp	SHORT $LN8@advanceAnd
$LN14@advanceAnd:
	mov	eax, DWORD PTR [rbx+72]
$LN8@advanceAnd:

; 419  : 
; 420  : 	// --- wrap if needed
; 421  : 	inputReadIndex &= wrapMask;

	mov	ecx, DWORD PTR [rbx+80]
	and	ecx, eax

; 422  : 
; 423  : 	return true;

	mov	al, 1
	mov	DWORD PTR [rbx+72], ecx
$LN1@advanceAnd:

; 424  : }

	add	rsp, 32					; 00000020H
	pop	rbx
	ret	0
?advanceAndCheckFFT@PhaseVocoder@@QEAA_NXZ ENDP		; PhaseVocoder::advanceAndCheckFFT
_TEXT	ENDS
; Function compile flags: /Ogspy
; File E:\prova\Shimmer\include\src\PSMVocoder.cpp
;	COMDAT ?processAudioSample@PhaseVocoder@@QEAANNAEA_N@Z
_TEXT	SEGMENT
this$ = 80
input$ = 88
fftReady$ = 96
?processAudioSample@PhaseVocoder@@QEAANNAEA_N@Z PROC	; PhaseVocoder::processAudioSample, COMDAT

; 437  : {

$LN8:
	mov	QWORD PTR [rsp+8], rbx
	push	rdi
	sub	rsp, 64					; 00000040H

; 438  : 	// --- if user did not manually do fft and overlap, do them here
; 439  : 	//     this allows maximum flexibility in use of the object
; 440  : 	if (needInverseFFT)

	cmp	BYTE PTR [rcx+96], 0
	mov	rdi, r8
	movaps	XMMWORD PTR [rsp+48], xmm6
	mov	rbx, rcx
	movaps	XMMWORD PTR [rsp+32], xmm7
	movaps	xmm7, xmm1
	je	SHORT $LN2@processAud

; 479  : 	fftw_execute(plan_backward);

	mov	rcx, QWORD PTR [rcx+32]
	call	QWORD PTR __imp_fftw_execute

; 480  : 
; 481  : 	// --- output is now in ifft_result array
; 482  : 	needInverseFFT = false;

	mov	BYTE PTR [rbx+96], 0
$LN2@processAud:

; 441  : 		doInverseFFT();
; 442  : 	if (needOverlapAdd)

	cmp	BYTE PTR [rbx+97], 0
	je	SHORT $LN3@processAud

; 443  : 		doOverlapAdd();

	xor	r8d, r8d
	xor	edx, edx
	mov	rcx, rbx
	call	?doOverlapAdd@PhaseVocoder@@QEAAXPEANH@Z ; PhaseVocoder::doOverlapAdd
$LN3@processAud:

; 444  : 
; 445  : 	fftReady = false;
; 446  : 
; 447  : 	// --- get the current output sample first
; 448  : 	double currentOutput = outputBuffer[outputReadIndex];

	mov	ecx, DWORD PTR [rbx+76]
	mov	rax, QWORD PTR [rbx+56]
	mov	BYTE PTR [rdi], 0
	movsd	xmm6, QWORD PTR [rax+rcx*8]

; 449  : 
; 450  : 	// --- set the buffer to 0.0 in preparation for the next overlap/add process
; 451  : 	outputBuffer[outputReadIndex++] = 0.0;

	and	QWORD PTR [rax+rcx*8], 0
	mov	ecx, DWORD PTR [rbx+76]

; 452  : 
; 453  : 	// --- wrap
; 454  : 	outputReadIndex &= wrapMaskOut;

	mov	eax, DWORD PTR [rbx+84]
	inc	ecx
	and	eax, ecx

; 455  : 
; 456  : 	// --- push into buffer
; 457  : 	inputBuffer[inputWriteIndex++] = (double)input;

	mov	ecx, DWORD PTR [rbx+64]
	mov	DWORD PTR [rbx+76], eax
	mov	rax, QWORD PTR [rbx+48]
	movsd	QWORD PTR [rax+rcx*8], xmm7
	mov	ecx, DWORD PTR [rbx+64]

; 458  : 
; 459  : 	// --- wrap
; 460  : 	inputWriteIndex &= wrapMask;

	mov	eax, DWORD PTR [rbx+80]
	inc	ecx
	and	eax, ecx

; 461  : 
; 462  : 	// --- check the FFT
; 463  : 	fftReady = advanceAndCheckFFT();

	mov	rcx, rbx
	mov	DWORD PTR [rbx+64], eax
	call	?advanceAndCheckFFT@PhaseVocoder@@QEAA_NXZ ; PhaseVocoder::advanceAndCheckFFT

; 464  : 
; 465  : 	return currentOutput;
; 466  : }

	mov	rbx, QWORD PTR [rsp+80]
	movaps	xmm0, xmm6
	movaps	xmm6, XMMWORD PTR [rsp+48]
	movaps	xmm7, XMMWORD PTR [rsp+32]
	mov	BYTE PTR [rdi], al
	add	rsp, 64					; 00000040H
	pop	rdi
	ret	0
?processAudioSample@PhaseVocoder@@QEAANNAEA_N@Z ENDP	; PhaseVocoder::processAudioSample
_TEXT	ENDS
; Function compile flags: /Ogspy
; File E:\prova\Shimmer\include\src\PSMVocoder.cpp
;	COMDAT ?doInverseFFT@PhaseVocoder@@QEAAXXZ
_TEXT	SEGMENT
this$ = 48
?doInverseFFT@PhaseVocoder@@QEAAXXZ PROC		; PhaseVocoder::doInverseFFT, COMDAT

; 477  : {

$LN4:
	push	rbx
	sub	rsp, 32					; 00000020H
	mov	rbx, rcx

; 478  : 	// do the IFFT
; 479  : 	fftw_execute(plan_backward);

	mov	rcx, QWORD PTR [rcx+32]
	call	QWORD PTR __imp_fftw_execute

; 480  : 
; 481  : 	// --- output is now in ifft_result array
; 482  : 	needInverseFFT = false;

	mov	BYTE PTR [rbx+96], 0

; 483  : }

	add	rsp, 32					; 00000020H
	pop	rbx
	ret	0
?doInverseFFT@PhaseVocoder@@QEAAXXZ ENDP		; PhaseVocoder::doInverseFFT
_TEXT	ENDS
; Function compile flags: /Ogspy
; File E:\prova\Shimmer\include\src\PSMVocoder.cpp
;	COMDAT ?doOverlapAdd@PhaseVocoder@@QEAAXPEANH@Z
_TEXT	SEGMENT
this$ = 8
outputData$ = 16
length$ = 24
?doOverlapAdd@PhaseVocoder@@QEAAXPEANH@Z PROC		; PhaseVocoder::doOverlapAdd, COMDAT

; 497  : 	// --- overlap/add with output buffer
; 498  : 	//     NOTE: this assumes input and output hop sizes are the same!
; 499  : 	outputWriteIndex = outputReadIndex;

	mov	eax, DWORD PTR [rcx+76]
	mov	r9, rcx
	mov	DWORD PTR [rcx+68], eax

; 500  : 
; 501  : 	if (outputData)

	test	rdx, rdx
	je	SHORT $LN8@doOverlapA

; 502  : 	{
; 503  : 		for (int i = 0; i < length; i++)

	movsxd	r10, r8d
	test	r8d, r8d
	jle	SHORT $LN6@doOverlapA
	xor	r8d, r8d
$LL4@doOverlapA:

; 504  : 		{
; 505  : 			// --- if you need to window the data, do so prior to this function call
; 506  : 			outputBuffer[outputWriteIndex++] += outputData[i];

	movsd	xmm0, QWORD PTR [rdx+r8*8]
	inc	r8
	mov	ecx, DWORD PTR [r9+68]
	mov	rax, QWORD PTR [r9+56]
	addsd	xmm0, QWORD PTR [rax+rcx*8]
	movsd	QWORD PTR [rax+rcx*8], xmm0
	mov	eax, DWORD PTR [r9+68]
	inc	eax

; 507  : 
; 508  : 			// --- wrap if index > bufferlength - 1
; 509  : 			outputWriteIndex &= wrapMaskOut;

	and	eax, DWORD PTR [r9+84]
	mov	DWORD PTR [r9+68], eax
	cmp	r8, r10
	jl	SHORT $LL4@doOverlapA

; 510  : 		}
; 511  : 		needOverlapAdd = false;
; 512  : 		return;

	jmp	SHORT $LN6@doOverlapA
$LN8@doOverlapA:

; 513  : 	}
; 514  : 
; 515  : 	for (int i = 0; i < frameLength; i++)

	xor	r8d, r8d
	cmp	DWORD PTR [rcx+104], r8d
	jbe	SHORT $LN6@doOverlapA
	xor	r10d, r10d
$LL7@doOverlapA:

; 516  : 	{
; 517  : 		// --- accumulate
; 518  : 		outputBuffer[outputWriteIndex++] += windowHopCorrection * ifft_result[i][0];

	mov	rax, QWORD PTR [r9+16]
	inc	r8d
	mov	edx, DWORD PTR [r9+68]
	mov	rcx, QWORD PTR [r9+56]
	movsd	xmm0, QWORD PTR [rax+r10]
	add	r10, 16
	mulsd	xmm0, QWORD PTR [r9+88]
	addsd	xmm0, QWORD PTR [rcx+rdx*8]
	movsd	QWORD PTR [rcx+rdx*8], xmm0
	mov	eax, DWORD PTR [r9+68]
	inc	eax

; 519  : 
; 520  : 		// --- wrap if index > bufferlength - 1
; 521  : 		outputWriteIndex &= wrapMaskOut;

	and	eax, DWORD PTR [r9+84]
	mov	DWORD PTR [r9+68], eax
	cmp	r8d, DWORD PTR [r9+104]
	jb	SHORT $LL7@doOverlapA
$LN6@doOverlapA:

; 522  : 	}
; 523  : 
; 524  : 	// --- set a flag
; 525  : 	needOverlapAdd = false;
; 526  : }

	mov	BYTE PTR [r9+97], 0
	ret	0
?doOverlapAdd@PhaseVocoder@@QEAAXPEANH@Z ENDP		; PhaseVocoder::doOverlapAdd
_TEXT	ENDS
END
