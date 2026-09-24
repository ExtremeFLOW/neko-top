!> @file nekotop_logger.f90
!! @copyright
!! Copyright (c) 2024-2026, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.

!> @brief Neko-TOP logging.
!! @details
!! This module declares the global Neko-TOP log object, reusing Neko's
!! `log_t` type. The log levels of Neko are re-exported here, such that a
!! `use nekotop_logger` alone provides both the log object and the levels.
module nekotop_logger
  use logger, only: log_t, LOG_SIZE, SEC_HEAD_SIZE, NEKO_LOG_QUIET, &
       NEKO_LOG_INFO, NEKO_LOG_VERBOSE, NEKO_LOG_DEPRECATION, NEKO_LOG_DEBUG
  implicit none
  private

  public :: log_t, LOG_SIZE, SEC_HEAD_SIZE, NEKO_LOG_QUIET, NEKO_LOG_INFO, &
       NEKO_LOG_VERBOSE, NEKO_LOG_DEPRECATION, NEKO_LOG_DEBUG

  !------------------------------------------------------------------
  ! Global log object
  !------------------------------------------------------------------

  !> Global Neko-TOP log stream.
  !! @details
  !! This is Neko-TOP's own log stream, independent of Neko's `neko_log`.
  !! It is initialized with `call nekotop_log%init(env_prefix = "NEKO_TOP")`,
  !! which reads `NEKO_TOP_LOG_LEVEL`, `NEKO_TOP_LOG_FILE` and
  !! `NEKO_TOP_LOG_TAB_SIZE`. These are separate from Neko's `NEKO_LOG_*`
  !! variables, so the verbosity and output destination of Neko-TOP can be
  !! controlled independently of Neko. It is freed with
  !! `call nekotop_log%free()`.
  type(log_t), public, target :: nekotop_log

end module nekotop_logger
