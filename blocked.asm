.section .text.spleen::blocked_expf64_strided_iter,"ax",@progbits
	.globl	spleen::blocked_expf64_strided_iter
	.p2align	4
.type	spleen::blocked_expf64_strided_iter,@function
spleen::blocked_expf64_strided_iter:
	.cfi_startproc
	push rbp
	.cfi_def_cfa_offset 16
	push r15
	.cfi_def_cfa_offset 24
	push r14
	.cfi_def_cfa_offset 32
	push r13
	.cfi_def_cfa_offset 40
	push r12
	.cfi_def_cfa_offset 48
	push rbx
	.cfi_def_cfa_offset 56
	sub rsp, 184
	.cfi_def_cfa_offset 240
	.cfi_offset rbx, -56
	.cfi_offset r12, -48
	.cfi_offset r13, -40
	.cfi_offset r14, -32
	.cfi_offset r15, -24
	.cfi_offset rbp, -16
	cmp rdx, rsi
	ja .LBB1_61
	mov rbx, rdx
	mov rax, rsi
	or rax, rdx
	shr rax, 32
	je .LBB1_2
	mov rax, rsi
	xor edx, edx
	div rbx
	jmp .LBB1_4
.LBB1_2:
	mov eax, esi
	xor edx, edx
	div ebx
.LBB1_4:
	mul rbx
	seto dl
	cmp rbx, 8
	setb r9b
	or r9b, dl
	cmp rax, rsi
	seta dl
	or dl, r9b
	mov r10b, 5
	je .LBB1_5
.LBB1_11:
	mov eax, r10d
	add rsp, 184
	.cfi_def_cfa_offset 56
	pop rbx
	.cfi_def_cfa_offset 48
	pop r12
	.cfi_def_cfa_offset 40
	pop r13
	.cfi_def_cfa_offset 32
	pop r14
	.cfi_def_cfa_offset 24
	pop r15
	.cfi_def_cfa_offset 16
	pop rbp
	.cfi_def_cfa_offset 8
	vzeroupper
	ret
.LBB1_5:
	.cfi_def_cfa_offset 240
	mov r10b, 1
	test rax, rax
	je .LBB1_11
	mov r15, rax
	or r15, rbx
	mov rsi, rax
	shr r15, 32
	je .LBB1_7
	xor edx, edx
	div rbx
	mov r9, rax
	cmp rcx, r9
	jb .LBB1_12
	jmp .LBB1_10
.LBB1_7:
	xor edx, edx
	div ebx
	mov r9d, eax
	cmp rcx, r9
	jae .LBB1_10
.LBB1_12:
	cmp rsi, 8
	jae .LBB1_14
	mov r10b, 1
	jmp .LBB1_11
.LBB1_14:
	mov rax, rsi
	vmovupd ymm3, ymmword ptr [rdi]
	vmovupd ymm4, ymmword ptr [rdi + 32]
	movzx ebp, r8b
	cmp bpl, 1
	mov r10b, 1
	je .LBB1_18
	cmp ebp, 2
	jne .LBB1_29
	mov r12, rdx
	mov rsi, rax
	test r15, r15
	je .LBB1_17
	xor edx, edx
	div rbx
	jmp .LBB1_23
.LBB1_18:
	test rcx, rcx
	je .LBB1_37
	lea r11, [rdi + 32]
	lea rsi, [8*rbx]
	inc r9
	imul r9, rbx
	vmovsd xmm2, qword ptr [rip + .LCPI1_0]
	mov r14, rcx
.LBB1_20:
	sub r9, rbx
	cmp rbx, r9
	ja .LBB1_37
	vmulsd xmm2, xmm0, xmm2
	vbroadcastsd ymm1, xmm2
	vfmadd231pd ymm3, ymm1, ymmword ptr [r11 - 32]
	vfmadd231pd ymm4, ymm1, ymmword ptr [r11]
	add r11, rsi
	dec r14
	jne .LBB1_20
	jmp .LBB1_37
.LBB1_29:
	test rcx, rcx
	je .LBB1_37
	lea rsi, [r9 - 1]
	imul rsi, rbx
	lea rsi, [rdi + 8*rsi]
	add rsi, 32
	inc r9
	imul r9, rbx
	vmovsd xmm2, qword ptr [rip + .LCPI1_0]
	xor r11d, r11d
	mov r14, rcx
	xor r12d, r12d
.LBB1_31:
	sub r12, rbx
	lea r13, [r9 + r12]
	cmp r13, rbx
	jb .LBB1_37
	vmulsd xmm2, xmm0, xmm2
	vbroadcastsd ymm1, xmm2
	vfmadd231pd ymm3, ymm1, ymmword ptr [rsi + 8*r11 - 32]
	vfmadd231pd ymm4, ymm1, ymmword ptr [rsi + 8*r11]
	mov r11, r12
	dec r14
	jne .LBB1_31
	jmp .LBB1_37
.LBB1_17:
	xor edx, edx
	div ebx
.LBB1_23:
	mov rax, rsi
	sub rsi, rdx
	lea r9, [rdi + 8*rsi]
	mov qword ptr [rsp + 24], rdi
	mov qword ptr [rsp + 32], rsi
	mov qword ptr [rsp + 40], r9
	mov qword ptr [rsp + 48], rdx
	mov qword ptr [rsp + 56], rbx
	mov qword ptr [rsp + 64], 1
	test rcx, rcx
	mov r10b, 1
	mov rdx, r12
	je .LBB1_37
	vmovsd xmm2, qword ptr [rip + .LCPI1_0]
	mov r11, rcx
.LBB1_25:
	dec r11
	mov qword ptr [rsp + 72], r11
	mov rsi, qword ptr [rsp + 64]
	test rsi, rsi
	jne .LBB1_33
	mov r13, qword ptr [rsp + 32]
	mov rsi, qword ptr [rsp + 56]
	sub r13, rsi
	jb .LBB1_37
	mov r9, qword ptr [rsp + 24]
	lea r14, [r9 + 8*rsi]
	mov qword ptr [rsp + 24], r14
	mov qword ptr [rsp + 32], r13
	cmp rsi, 7
	jbe .LBB1_28
.LBB1_36:
	vmulsd xmm2, xmm0, xmm2
	vbroadcastsd ymm1, xmm2
	vfmadd231pd ymm3, ymm1, ymmword ptr [r9]
	vfmadd231pd ymm4, ymm1, ymmword ptr [r9 + 32]
	test r11, r11
	jne .LBB1_25
	jmp .LBB1_37
.LBB1_33:
	vmovapd xmmword ptr [rsp + 96], xmm2
	vmovupd ymmword ptr [rsp + 112], ymm4
	vmovupd ymmword ptr [rsp + 144], ymm3
	vmovaps xmmword ptr [rsp + 80], xmm0
	mov dword ptr [rsp + 4], r8d
	mov qword ptr [rsp + 16], rcx
	mov qword ptr [rsp + 8], rax
	mov r13, rdi
	mov qword ptr [rsp + 64], 0
	lea rdi, [rsp + 24]
	vzeroupper
	call <core::iter::adapters::map::Map<core::slice::iter::ChunksExact<f64>, <spleen::iter_utils::BlockedStridedSliceMut<8, f64> as spleen::iter_utils::StridedIter>::strided_iter::{closure#0}> as core::iter::traits::iterator::Iterator>::nth
	test rax, rax
	je .LBB1_34
	mov r9, rax
	mov r11, qword ptr [rsp + 72]
	mov rdi, r13
	mov rax, qword ptr [rsp + 8]
	mov rcx, qword ptr [rsp + 16]
	mov r8d, dword ptr [rsp + 4]
	mov r10b, 1
	mov rdx, r12
	vmovapd xmm0, xmmword ptr [rsp + 80]
	vmovupd ymm3, ymmword ptr [rsp + 144]
	vmovupd ymm4, ymmword ptr [rsp + 112]
	vmovapd xmm2, xmmword ptr [rsp + 96]
	jmp .LBB1_36
.LBB1_61:
	lea rdi, [rip + .Lanon.430998381b82c322518bf48d687a6884.4]
	call qword ptr [rip + core::option::unwrap_failed@GOTPCREL]
.LBB1_28:
	lea rdi, [rip + .Lanon.430998381b82c322518bf48d687a6884.1]
	vzeroupper
	call qword ptr [rip + core::option::unwrap_failed@GOTPCREL]
.LBB1_34:
	mov rdi, r13
	mov rax, qword ptr [rsp + 8]
	mov rcx, qword ptr [rsp + 16]
	mov r8d, dword ptr [rsp + 4]
	mov r10b, 1
	mov rdx, r12
	vmovaps xmm0, xmmword ptr [rsp + 80]
	vmovups ymm3, ymmword ptr [rsp + 144]
	vmovups ymm4, ymmword ptr [rsp + 112]
.LBB1_37:
	vmovupd ymmword ptr [rdi], ymm3
	vmovupd ymmword ptr [rdi + 32], ymm4
	mov rsi, rax
	sub rsi, rdx
	sub rsi, rbx
	jb .LBB1_11
	cmp rsi, rbx
	jb .LBB1_41
	vbroadcastsd ymm1, xmm0
	lea r9, [8*rbx]
	lea r10, [rdi + 8*rbx]
	add r10, 32
	mov r11, rsi
	.p2align	4
.LBB1_40:
	sub r11, rbx
	vfmadd213pd ymm3, ymm1, ymmword ptr [r10 - 32]
	vmovupd ymmword ptr [r10 - 32], ymm3
	vfmadd213pd ymm4, ymm1, ymmword ptr [r10]
	vmovupd ymmword ptr [r10], ymm4
	add r10, r9
	cmp rbx, r11
	jbe .LBB1_40
.LBB1_41:
	lea r9, [rdi + 8*rdx]
	vmovupd ymm1, ymmword ptr [r9 + 8*rsi]
	vmovupd ymm4, ymmword ptr [r9 + 8*rsi + 32]
	cmp r8b, 1
	je .LBB1_54
	cmp ebp, 2
	jne .LBB1_43
	mov rcx, rsi
	sub rcx, rbx
	jae .LBB1_57
.LBB1_10:
	xor r10d, r10d
	jmp .LBB1_11
.LBB1_54:
	vaddsd xmm2, xmm0, qword ptr [rip + .LCPI1_2]
	jmp .LBB1_55
.LBB1_43:
	test rcx, rcx
	je .LBB1_50
	mov r11, rax
	test r15, r15
	je .LBB1_45
	xor edx, edx
	div rbx
	mov rdx, rax
	jmp .LBB1_47
.LBB1_57:
	vbroadcastsd ymm2, xmm0
	vfmadd231pd ymm1, ymm2, ymmword ptr [r9 + 8*rcx]
	vfmadd231pd ymm4, ymm2, ymmword ptr [r9 + 8*rcx + 32]
	vmulsd xmm2, xmm0, xmm0
	vaddsd xmm2, xmm2, qword ptr [rip + .LCPI1_2]
.LBB1_55:
	vdivsd xmm2, xmm0, xmm2
	vbroadcastsd ymm2, xmm2
	vmulpd ymm1, ymm2, ymm1
	vmulpd ymm3, ymm2, ymm4
	jmp .LBB1_51
.LBB1_45:
	xor edx, edx
	div ebx
	mov edx, eax
.LBB1_47:
	inc rdx
	imul rdx, rbx
	lea r10, [rdi + 32]
	lea r8, [8*rbx]
	vmovsd xmm2, qword ptr [rip + .LCPI1_0]
	mov rax, r11
.LBB1_48:
	sub rdx, rbx
	cmp rbx, rdx
	ja .LBB1_50
	vmulsd xmm2, xmm0, xmm2
	vbroadcastsd ymm3, xmm2
	vfmadd231pd ymm1, ymm3, ymmword ptr [r10 - 32]
	vfmadd231pd ymm4, ymm3, ymmword ptr [r10]
	add r10, r8
	dec rcx
	jne .LBB1_48
.LBB1_50:
	vxorpd xmm2, xmm0, xmmword ptr [rip + .LCPI1_1]
	vbroadcastsd ymm2, xmm2
	vmulpd ymm1, ymm1, ymm2
	vmulpd ymm3, ymm4, ymm2
.LBB1_51:
	vmovupd ymmword ptr [r9 + 8*rsi], ymm1
	vmovupd ymmword ptr [r9 + 8*rsi + 32], ymm3
	mov r10b, -1
	cmp rsi, rbx
	jb .LBB1_11
	test r15, r15
	je .LBB1_53
	xor edx, edx
	div rbx
	jmp .LBB1_59
.LBB1_53:
	xor edx, edx
	div ebx
.LBB1_59:
	vbroadcastsd ymm2, xmm0
	add rax, -2
	imul rax, rbx
	lea rax, [rdi + 8*rax]
	add rax, 32
	xor ecx, ecx
	mov r10b, -1
	.p2align	4
.LBB1_60:
	vsubpd ymm1, ymm1, ymmword ptr [rax + 8*rcx - 32]
	vmulpd ymm1, ymm2, ymm1
	vmovupd ymmword ptr [rax + 8*rcx - 32], ymm1
	vsubpd ymm0, ymm3, ymmword ptr [rax + 8*rcx]
	vmulpd ymm3, ymm2, ymm0
	vmovupd ymmword ptr [rax + 8*rcx], ymm3
	sub rcx, rbx
	lea rdx, [rsi + rcx]
	cmp rdx, rbx
	jae .LBB1_60
	jmp .LBB1_11
