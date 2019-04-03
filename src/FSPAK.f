c      ----------------
       subroutine fspak
c      ----------------
     +                  (option,n,ia,ja,a,b,flag,iout,ioor,
     +                   available,needed,is,fnode,bnode,
     +                   fnseg,bnseg,fx,feqb,irank)
c -----------------------------------------------------------
c Public Sparse Matrix Package for Symmetric Matrices
c Authors: Miguel Perez-Enciso, Ignacy Misztal, Mauricio Elzo
c Mon Aug 24, 1992 12:39:46 - Jan 17, 1994
c -----------------------------------------------------------
      implicit none
      integer zero,order,symfac,numfac,solve,sinit,ssolve
     +       ,det,ldet,check,spin,spin1,restart,stat,equal,mx_st
      parameter (order=10,  symfac=20,  numfac=40, solve=50,
     +           sinit=51,  ssolve=52,  det=54,    ldet =55,
     +           check=56,  spin=60,  spin1=61, restart=70,
     +           stat=80,  equal=1,   mx_st=50 )
      integer available,r_avail
     +       ,delta,flag,ioor,iout,n,ia(1),ja(1),is(available)
     +       ,fnode(1),bnode(1),fnseg,bnseg,feqb,irank,ratio,fratio
     +       ,path,found,fpath,fstart,fpnode,bpath,bstart,bpnode
     +       ,maxju,maxu,needed,maxneed,i,adjncy,d,invp,iu,iju,ju,u
     +       ,maxint,nza,nzad,option,perm,iap,jap,ap,xadj,first,ios
     +       ,zja,za,zd,utia,utja,work1,work2,work3,work4,work5,work6
     +       ,even,count_na
      real*8 a(1),b(1),fx(1),st(mx_st),t_t,t0,tol,second
      equivalence (st(1),maxju),(st(2),maxu),(st(3),maxneed),(st(4),nza)
c     this is to avoid that the value of ratio is lost
c     with dynamic allocation memory
      save ratio
      data st/mx_st*0./,first/0/
c..   tolerance for accuracy in solution (in options numfac and check)
     +    ,tol/1.d-9/


c function to ensure that real*8 variables are assigned an even address.
      even(i)=2*(i/2)+1
      
      flag=0
c      t_t=second()
c     -----------------------------------------------
c     first time it checks whether ratio equal 1 or 2
c     -----------------------------------------------
      if (first.eq.0) then
         ratio=fratio()
         first=1
      endif
c     -------------------------------------------
c     no. of off-diagonal nze in upper triangular
c     -------------------------------------------
c ----comment either line; first one is faster but takes more memory
      nza=(ia(n+1)-1)
c      nza=count_nza(n,ia,ja,iout)
c     --------------------
c     no. of nze in adjncy
c     --------------------
      nzad=2*nza
c     ----------------
c     allocate storage
c     ----------------
      perm   = 1
      invp   = perm + n
c     -------------------------
      if (option.eq.order) then
c     -------------------------
         xadj   = invp   + n
         adjncy = xadj   + n+1
         work1  = adjncy + nzad
         work2  = work1 + n
         work3  = work2 + n
         work4  = work3 + n
         needed = work4 + n
         call checkst (option,needed,maxneed,available,iout,flag)
         if (flag.ne.0) return
c        --------------------------------------------------
c        max integer (used as flag in the ordering routine)
c        --------------------------------------------------
         maxint=9999999
c        ---------------------------------------------------------
c        parameter for type of ordering -1=MD, 0=MMD, >0=slack MMD
c        ---------------------------------------------------------
         delta=0
         call unfold
     +               (n,ia,ja,is(xadj),is(adjncy),is(work1))
c         t0=second()
c        -------------------------------
c        if mmd ordering routines by Liu
c        -------------------------------
         call genmmd
     +             (n,is(xadj),is(adjncy),is(invp),is(perm),delta,
     +              is(work1),is(work2),is(work3),is(work4),maxint,
     +              maxju)
c        ---------------------------
c        elseif md sparspak routines
c        ---------------------------
c        work5 = work4 + n
c        work6 = work5 + n
c        needed= work6 + n
c        call checkst (option,needed,maxneed,available,iout,flag)
c        if(flag.ne.0) return
c        call genqmd
c    +               (n,is(xadj),is(adjncy),is(perm),is(invp),
c    +                is(work1),is(work2),is(work3),is(work4),
c    +                is(work5),is(work6),maxju)
c         st(21)=st(21)+second()-t0
c        ------------------------------
c        write down permutation vectors
c        ------------------------------
         rewind ioor
         write(ioor,*) n,maxju
         write(ioor,'(10i7)') (is(i),i=perm,perm+n-1),
     +                        (is(i),i=invp,invp+n-1)
c     -------------------------------
      elseif (option.eq.restart) then
c     -------------------------------
         rewind ioor
         i=0
         read(ioor,*,end=1) i,maxju
1        if(i.eq.n) then
            read(ioor,'(10i7)',end=2) (is(i),i=perm,perm+n-1),
     +                                (is(i),i=invp,invp+n-1)
2           inquire(ioor,iostat=ios)
            if (ios.ne.0) then
               write(iout,*) 'error in unit IOOR '
               flag=-2
            endif
         else
            if (i.ne.0) then
            write(iout,*) 'incorrect size of permutation vector'
            write(iout,*) 'size ',i,' found; size ',n, ' required'
            endif
            flag=-1
         endif
c     ------------------------------
      elseif (option.eq.symfac) then
c     ------------------------------
         iu     = invp   + n
         iju    = iu     + n+1
         ju     = iju    + n+1
         xadj   = ju     + maxju
         adjncy = xadj   + n+1
         work1  = adjncy + nzad
         work2  = work1  + n
         work3  = work2  + n
         needed = work3  + n
         call checkst (option,needed,maxneed,available,iout,flag)
         if (flag.ne.0) return
c         t0=second()
         call unfold
     +               (n,ia,ja,is(xadj),is(adjncy),is(work1))
         call smbfct
     +               (n,is(xadj),is(adjncy),is(perm),is(invp),is(iu),
     +                maxu,is(iju),is(ju),maxju,is(work1),is(work2),
     +                is(work3),flag)
c         st(22)=st(22)+second()-t0
      else
         iu        = invp  + n
         iju       = iu    + n+1
         ju        = iju   + n+1
         d         = even(ju    + maxju)
         u         = even(d     + n*ratio)
c        --------------------------
         if (option.eq.numfac) then
c        --------------------------
            iap    = u     + maxu*ratio
            jap    = iap   + n+1
            ap     = even(jap   + nza)
            work1  = even(ap    + nza*ratio)
            needed = work1 + n*ratio
            call checkst (option,needed,maxneed,available,iout,flag)
            if (flag.ne.0) return
c            t0=second()
c           ---------------------------------------
c           first load numerical values of a into u
c           ---------------------------------------
            call loadap
     +                (n,ia,ja,a,is(iap),is(jap),is(ap),is(d),is(invp),
     +                 is(work1))
            call loadu
     +                (n,is(iap),is(jap),is(ap),is(iu),is(ju),is(iju),
     +                 is(u),is(work1))
                             
            work1  = u     + maxu*ratio
            work2  = work1 + n
            work3  = even(work2 + n)
            needed = work3 + n*ratio
            call checkst (option,needed,maxneed,available,iout,flag)
            if (flag.ne.0) return
c           -----------------------
c           numerical factorization (slightly modified sparspak routine)
c           ----------------------------------
            call fgsfct
     +                 (n,is(iu),is(u),is(iju),is(ju),is(d),is(work1),
     +                  is(work2),is(work3),flag,tol,irank)
            call merlin (n,is(iu),is(u),is(d),tol)
c            st(24)=st(24)+second()-t0
            st(11)=st(11)+1
c           ----------
c           zero pivot
c           ----------
            if(flag.ne.0) flag=numfac*1000000 + flag
c        -----------------------------
         elseif (option.eq.solve) then
c        -----------------------------
            work1  = even(u     + maxu*ratio)
            needed = work1 + n*ratio
c            t0=second()

            call permute (n,b,is(work1),is(invp))
            call fullfb
     +                   (n,is(iu),is(ju),is(iju),is(u),is(d),is(work1))
     
            call permute (n,is(work1),b,is(perm))
c            st(25)=st(25)+second()-t0
            st(12)=st(12)+1
c        --------------------------------------------------
         elseif (option.ge.sinit.and.option.le.ssolve) then
c        --------------------------------------------------
            path  = u     + maxu*ratio
            found = path  + n
            work1 = even(found + n)
            fstart= work1 + n*ratio
            fpnode= fstart+ fnseg+1
            needed= fpnode+ fnseg
            call checkst
c                        ------------------------------------------
c                        allows a safety margin for fpath and bpath
c                        ------------------------------------------
     +                   (option,needed+200,maxneed,available,iout,
     +                    flag)
            if(flag.ne.0) return
c           -------------------------
            if (option.eq.sinit) then
c           -------------------------
c              initialize and computes factorization table
c              -------------------------------------------
               call rzero (n,is(work1))
               call lzero (n,is(found))
               call makepath (n,is(path),is(iu),is(ju),is(iju))
c           ------------------------------
            elseif (option.eq.ssolve) then
c           ------------------------------
c               t0=second()
               fpath = fpnode + fnseg
c              ------------------------------------------------------
c              obtains forward path (in the permuted vector of nodes)
c              ------------------------------------------------------
               call spermute (fnseg,fnode,is(fpnode),is(invp))
               call segpath
     +                      (n,is(fpnode),fnseg,is(fstart),is(fpath),
     +                       is(path),is(found))
               st(34)=st(34)+is(fstart+fnseg)-1
               if(feqb.eq.equal) then
                  needed= fpath + is(fstart+fnseg)-1
                  call checkst
     +                        (option,needed,maxneed,available,iout,
     +                         flag)
                  if (flag.ne.0) return
c                 -----------------------------
c                 initialize working rhs vector
c                 -----------------------------
                  call zerout (is(work1),is(fpath),is(fstart),fnseg)
c                 ----------------
c                 input rhs values
c                 ----------------
                  call expand (fnseg,is(fpnode),fx,is(work1))
c                 -----
c                 solve
c                 -----
                  call fastfb
     +                        (is(fpath),is(fstart),fnseg,is(fpath),
     +                         is(fstart),fnseg,is(iu),is(ju),is(iju),
     +                         is(u),is(d),is(work1))
c                 ---------------------
c                 stores solutions in b
c                 ---------------------
                  call compress (fnseg,is(fpnode),b,is(work1))
                  st(31)= max0(int(st(31)),is(fstart+fnseg)-1)
                  st(35)=st(35)+is(fstart+fnseg)-1
               else
                  bstart = fpath + is(fstart+fnseg)-1
                  bpnode = bstart+ bnseg+1
                  bpath  = bpnode+ bnseg
c                 ----------------------------------
c                 obtains backward path if necessary
c                 ----------------------------------
                  call spermute (bnseg,bnode,is(bpnode),is(invp))
                  call segpath
     +                         (n,is(bpnode),bnseg,is(bstart),
     +                          is(bpath),is(path),is(found))
                  needed= bpath + is(bstart+bnseg)-1
                  call checkst
     +                        (option,needed,maxneed,available,iout,
     +                         flag)
                  if (flag.ne.0) return
                  call zerout (is(work1),is(bpath),is(bstart),bnseg)
                  call zerout (is(work1),is(fpath),is(fstart),fnseg)
                  call expand (fnseg,is(fpnode),fx,is(work1))
                  call fastfb
     +                        (is(fpath),is(fstart),fnseg,is(bpath),
     +                         is(bstart),bnseg,is(iu),is(ju),is(iju),
     +                         is(u),is(d),is(work1))
                  call compress (bnseg,is(bpnode),b,is(work1))
                  st(31)=max(int(st(31)),is(bstart+bnseg)-1)
                  st(35)=st(35)+is(bstart+bnseg)-1
               endif
               st(13)=st(13)+1
c               st(27)=st(27)+second()-t0
               st(30)=max(int(st(30)),fnseg)
               st(31)=max(int(st(31)),bnseg)
               st(32)=st(32)+fnseg
               st(33)=st(33)+bnseg
            endif
c        ---------------------------
         elseif (option.eq.det) then
c        ---------------------------
            call mkdet (n,is(d),b,irank,tol,iout,flag)
            st(14)=st(14)+1
c        ----------------------------
         elseif (option.eq.ldet) then
c        ----------------------------
            call mkldet (n,is(d),b,irank,tol,iout,flag)
            st(14)=st(14)+1
c        ----------------------------
         elseif (option.eq.check) then
c        ----------------------------
            work1  = u     + maxu*ratio
            call checkb 
     +                  (n,ia,ja,a,is(d),b,fx,tol,is(work1),is(invp)
     +                  ,iout,flag)
c        ----------------------------
         elseif (option.eq.spin) then
c        ----------------------------
            zja    = u     + maxu*ratio
            za     = even(zja   + maxu+n)
            zd     = even(za    + (maxu+n)*ratio)
            utia   = zd    + n*ratio
            utja   = utia  + n+1
            work1  = utja  + maxu+n
            work2  = work1 + n
            work3  = even(work2 + n)
            needed = work3 + n*ratio
            call checkst (option,needed,maxneed,available,iout,flag)
            if (flag.ne.0) return
c            t0=second()
            call fsinverse
     +                    (n,is(perm),is(invp),is(d),is(iju),is(ju),
     +                     is(iu),is(u),is(work3),is(za),is(zd),
     +                     is(zja),is(utia),is(utja),
     +                     is(work1),is(work2),ia,ja,a)
c            st(26)=st(26)+second()-t0
            st(15)=st(15)+1
c        ----------------------------
         elseif (option.eq.spin1) then
c        ----------------------------
            work1  = even(u     + maxu*ratio)
            work2  = even(work1 + n*ratio)
            work5  = even(work2 + n*ratio)
            needed = work5 + nza+n
            call checkst (option,needed,maxneed,available,iout,flag)
            if (flag.ne.0) return
c            t0=second()
            call fsinverse1
     +                    (n,is(perm),is(invp),is(d),is(iju),is(ju),
     +                     is(iu),is(u),is(work1),is(work2),is(work1),
     +                     is(work2),is(work5),ia,ja,a)
c            st(26)=st(26)+second()-t0
            st(15)=st(15)+1
c        ----------------------------
         elseif (option.eq.stat) then
c        ----------------------------
            call figureout (n,available,iout,st,irank)
c        ----
         else
c        ----
            write(iout,*) 'option no. ',option,' not implemented'
            flag=-99
         endif
      endif
c      st(20)=st(20)+second()-t_t
      return
      end


c last modified Nov 23, 1994
c
c list of subroutines needed for FSPAK
c Miguel PEREZ-ENCISO
c Ignacy MISZTAL
c Mauricio ELZO
c
c second
c count_nza
c checkb
c checkst
c fastfb
c fullfb
c fgsfct
c figureout
c fratio
c fsinverse
c fsinverse1
c ival
c loadap
c loadu
c lzero
c makepath
c mkdet
c mkldet
c permute
c rowtocol
c rowtocolu
c rzero
c segpath
c spermute
c unfold
c zero
c zerout
c
      function count_nza(n,ia,ja,iout)
c     mpe, 17 January 1995
c counts no. of nze off-diagonal and checks whether there
c is any empty row.
      implicit none
      integer count_nza,n,i,iout,ia(1),ja(1),k,n0
      n0=0
      count_nza=0
      do i=1,n
         if(ia(i).eq.ia(i+1)) then
            n0=n0+1
         else
            do k=ia(i),ia(i+1)-1
               if(ja(k).ne.i) count_nza=count_nza+1
            enddo
         endif
      enddo
      if (n0.ne.0) write(iout,*) 'WARNING! ',n0,' empty rows'
      return
      end
c     -----------------
      subroutine checkb (n,ia,ja,a,d,b,x,tol,x1,perm,iout,flag)
c     -----------------
c     mpe, 18 Dec 1992 - 17 January 1995
      implicit none
      integer n,ia(1),ja(1),iout,flag,i,j,k,perm(1)
      real*8 dev,tol,a(1),b(1),d(1),x(1),x1(1)
      dev=0.
      do i=1,n
         x1(i)=0.
      enddo
      do i=1,n
         do k=ia(i),ia(i+1)-1
            j=ja(k)
            x1(i)=x1(i)+a(k)*b(j)
            if(i.lt.j) x1(j)=x1(j)+a(k)*b(i)
         enddo
      enddo
      do i=1,n
c         WRITE(IOUT,*) I,X1(I),X(I)
         if(ia(i+1).gt.ia(i)) dev=dev+abs(x1(i)-x(i))/n
      enddo
      if(dev.gt.tol) then
          flag=1
          write(iout,*) 'ERROR in CHECKB: not accurate solution'
          write(iout,*) 'observed deviation = ',dev
          write(iout,*) 'maximum allowed    = ',tol
C CHANGED 11/25/96          do i=1,n
C            write(iout,*) i,x1(i),x(i)
C          end do
      endif
c     kgb - return dev in position 1 of x(1)=fx(1)
      x(1)=dev
      return
      end

c     ------------------
      subroutine checkst (option,needed,maxneed,available,iout,flag)
c     ------------------
      integer option,needed,maxneed,available,iout,flag
      maxneed=max(maxneed,needed)
      if(maxneed.gt.available) then
         write(iout,*) ' insufficient storage in option ',option
         write(iout,*) ' avaliable ',available,'; needed ',needed
         flag=option*100+1
      endif
      return
      end

c     -------------------
      subroutine compress (nc,ipos,bc,b)
c     -------------------
c compreses nc positions of b into bc
      integer ipos(1)
      real*8 bc(1),b(1)
      do i=1,nc
         bc(i)=b(ipos(i))
      enddo
      return
      end

c     -----------------
      subroutine expand (nc,ipos,bc,b)
c     -----------------
c expand a compressed vector bc into b
c MPE
      integer ipos(1)
      real*8 bc(1),b(1)
      do i=1,nc
         b(ipos(i))=bc(i)
      enddo
      return
      end

c     -----------------
      subroutine fastfb
c     -----------------
     +                  (fpath,fstart,fnseg,bpath,bstart,bnseg,
     +                   iu,ju,iju,u,d,b)
c performs fast forward/backward on a triangular system iuju,iju,u
c M_P-E Thu Aug 13, 1992 14:34:41
      implicit none
      integer fpath(1),fstart(1),fnseg,bpath(1),bstart(1),bnseg,
     +        iu(1),ju(1),iju(1),i,j,k,ik,iseg
      real*8 u(1),d(1),b(1),bi
c fast forward in order within segments with reverse order btw segments
      do iseg=fnseg,1,-1
         do ik=fstart(iseg),fstart(iseg+1)-1
            i=fpath(ik)
            bi=b(i)
            do k=iu(i),iu(i+1)-1
               j=ju(iju(i)+k-iu(i))
               b(j)=b(j)+u(k)*bi
            enddo
          enddo
c fast diagonal
          do ik=fstart(iseg),fstart(iseg+1)-1
             i=fpath(ik)
             b(i)=b(i)*d(i)
          enddo
       enddo
c fast backward in reverse order within segments with reverse order btw seg
c segm
       do iseg=1,bnseg
          do ik=bstart(iseg+1)-1,bstart(iseg),-1
             i=bpath(ik)
             bi=b(i)
             do k=iu(i),iu(i+1)-1
                j=ju(iju(i)+k-iu(i))
                bi=bi+u(k)*b(j)
             enddo
             b(i)=bi
          enddo
       enddo
       return
       end

c     -----------------
      subroutine fullfb (n,iu,ju,iju,u,d,b)
c     -----------------
c performs full forward/backward on a triangular system iuju,iju,u
c M_P-E Sat Jul 11, 1992 17:54:59
      implicit none
      integer n,iu(1),ju(1),iju(1),i,j,k
      real*8 u(1),d(1),b(1),bi
c full forward
      do i=1,n         
         bi=b(i)
         do k=iu(i),iu(i+1)-1
            j=ju(iju(i)+k-iu(i))
            b(j)=b(j)+u(k)*bi
         enddo
      enddo
      
c diagonal
       do i=1,n
          b(i)=b(i)*d(i)
       enddo
       
c full backward
       do i=n,1,-1
          bi=b(i)
          do k=iu(i),iu(i+1)-1
             j=ju(iju(i)+k-iu(i))
             bi=bi+u(k)*b(j)
          enddo
          b(i)=bi
       enddo
       return
       end

c     ---------------
      function fratio()
c     ---------------
c returns the ratio of storage occupied by reals and integers;
c the only two values returned are 1 and 2. Orginally written to
c facilitate easy porting between Crays (ratio=1) to other
c computers (ratio=2, if reals are real*8)
c Ignacy Misztal
      integer fratio,ix(2),const
      parameter(const=12345)
      real*8 y
      equivalence (y,ix(1))
      ix(2)=const
      y=0
      if (ix(2).eq.const) then
         fratio=1
        else
         fratio=2
      endif
      end

csppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
c Modified for sp semidefinite matrices (mae,1993)
c (partial inclusion of modifications done by Kachman (UNE)
c  on subr. gsfct as described by Curt Finley (UCD), 1993)
C Slightly modified and restructured
C M_P-E Sat Jul 11, 1992 16:50:38
C From George & Liu (1981)
C----- SUBROUTINE FGSFCT
C***************************************************************
C***************************************************************
C******     GSFCT ..... GENERAL SPARSE SYMMETRIC FACT     ******
C***************************************************************
C***************************************************************
C
C     PURPOSE - THIS SUBROUTINE PERFORMS THE SYMMETRIC
C        FACTORIZATION FOR A GENERAL SPARSE SYSTEM, STORED IN
C        THE COMPRESSED SUBSCRIPT DATA FORMAT.
C                                                                         1
C     INPUT PARAMETERS -                                                  1
C        NEQNS - NUMBER OF EQUATIONS.                                     1
C        XLNZ - INDEX VECTOR FOR LNZ.  XLNZ(I) POINTS TO THE              1
C               START OF NONZEROS IN COLUMN I OF FACTOR L.                1
C        (XNZSUB, NZSUB) - THE COMPRESSED SUBSCRIPT DATA                  1
C               STRUCTURE FOR FACTOR L.                                   1
C                                                                         1
C     UPDATED PARAMETERS -                                                1
C        LNZ - ON INPUT, CONTAINS NONZEROS OF A, AND ON                   1
C               RETURN, THE NONZEROS OF L.                                2
C        DIAG - THE DIAGONAL OF L OVERWRITES THAT OF A.                   2
C        IFLAG - THE ERROR FLAG.  IT IS SET TO 1 IF A ZERO OR             2
C               NEGATIVE SQUARE ROOT OCCURS DURING THE                    2
C               FACTORIZATION.                                            2
C        OPS   - A DOUBLE PRECISION COMMON PARAMETER THAT IS              2
C                INCREMENTED BY THE NUMBER OF OPERATIONS                  2
C                PERFORMED BY THE SUBROUTINE.                             2
C                                                                         2
C     WORKING PARAMETERS -                                                2
C        LINK - AT STEP J, THE LIST IN                                    3
C                  LINK(J), LINK(LINK(J)), ...........                    3
C               CONSISTS OF THOSE COLUMNS THAT WILL MODIFY                3
C               THE COLUMN L(*,J).                                        3
C        FIRST - TEMPORARY VECTOR TO POINT TO THE FIRST                   3
C               NONZERO IN EACH COLUMN THAT WILL BE USED                  3
C               NEXT FOR MODIFICATION.                                    3
C        TEMP - A TEMPORARY VECTOR TO ACCUMULATE MODIFICATIONS.           3
C                                                                         3
C***************************************************************          3
C                                                                         4
cmae
      SUBROUTINE  FGSFCT ( NEQNS, XLNZ, LNZ, XNZSUB, NZSUB, DIAG,
     1                     LINK, FIRST, TEMP, IFLAG,                      4
     2                     tol,irank)
C                                                                         4
C***************************************************************          4
C                                                                         4
cmae
         integer*4 irank
         real*8 tol
         REAL*8 COUNT, OPS
         COMMON  /SPKOPS/ OPS                                             4
         REAL*8 DIAG(1), LNZ(1), TEMP(1), DIAGJ, LJK
         INTEGER LINK(1), NZSUB(1)                                        4
         INTEGER FIRST(1), XLNZ(1), XNZSUB(1),                            5
     1           I, IFLAG, II, ISTOP, ISTRT, ISUB, J,                     5
     1           K, KFIRST, NEQNS, NEWK                                   5
C                                                                         5
C***************************************************************          5
C                                                                         5
cmae
c        ------------------------------
c        initialize irank and dmin
c        ------------------------------
         irank=0
         dmin=1.d-08
cmae
C        ------------------------------                                   5
C        INITIALIZE WORKING VECTORS ...                                   5
C        ------------------------------                                   5
         DO I = 1, NEQNS
            LINK(I) = 0                                                   6
            TEMP(I) = 0.0                                                 6
         ENDDO
C        --------------------------------------------                     6
C        COMPUTE COLUMN L(*,J) FOR J = 1,...., NEQNS.                     6
C        --------------------------------------------                     6
         DO J = 1, NEQNS
C           -------------------------------------------                   6
C           FOR EACH COLUMN L(*,K) THAT AFFECTS L(*,J).                   6
C           -------------------------------------------                   6
            DIAGJ = 0.0                                                   7
            NEWK = LINK(J)                                                7
            K    = NEWK
            DO WHILE (K.NE.0)
               NEWK = LINK(K)                                             7
C              ---------------------------------------                    7
C              OUTER PRODUCT MODIFICATION OF L(*,J) BY                    7
C              L(*,K) STARTING AT FIRST(K) OF L(*,K).                     7
C              ---------------------------------------                    7
               KFIRST = FIRST(K)                                          7
               LJK    = LNZ(KFIRST)                                       8
               DIAGJ = DIAGJ + LJK*LJK                                    8
               OPS  = OPS + 1.00                                          8
               ISTRT = KFIRST + 1                                         8
               ISTOP = XLNZ(K+1) - 1                                      8
               IF ( ISTOP .GE. ISTRT )  THEN
C                 ------------------------------------------              8
C                 BEFORE MODIFICATION, UPDATE VECTORS FIRST,              8
C                 AND LINK FOR FUTURE MODIFICATION STEPS.                 8
C                 ------------------------------------------              8
                  FIRST(K) = ISTRT                                        9
                  I = XNZSUB(K) + (KFIRST-XLNZ(K)) + 1                    9
                  ISUB = NZSUB(I)                                         9
                  LINK(K) = LINK(ISUB)                                    9
                  LINK(ISUB) = K                                          9
C                 ---------------------------------------                 9
C                 THE ACTUAL MOD IS SAVED IN VECTOR TEMP.                 9
C                 ---------------------------------------                 9
                  DO II = ISTRT, ISTOP                                    9
                     ISUB = NZSUB(I)                                      9
                     TEMP(ISUB) = TEMP(ISUB) + LNZ(II)*LJK               10
                     I = I + 1                                           10
                  ENDDO                                                  10
                  COUNT = ISTOP - ISTRT + 1                              10
                  OPS  = OPS + COUNT
               ENDIF
               K    = NEWK
            ENDDO
C           ----------------------------------------------               10
C           APPLY THE MODIFICATIONS ACCUMULATED IN TEMP TO               10
C           COLUMN L(*,J).                                               10
C           ----------------------------------------------               10
            DIAGJ = DIAG(J) - DIAGJ                                      11
cmae
c           IF ( DIAGJ .LE. 0.0E0 )  THEN
C              ------------------------------------------------------
C              ERROR - ZERO OR NEGATIVE SQUARE ROOT IN FACTORIZATION.
C              ------------------------------------------------------
c              IFLAG = J
c              RETURN
c           ENDIF
c           DIAGJ = SQRT(DIAGJ)
cmae        DIAG(J) = DIAGJ
cmae
c  this is the core of Kachman's modification to subr. gsfct
            if (diag(j).ge.dmin) then
               if (diagj.ge.(tol*diag(j))) then
                  irank=irank+1
                  diagj = sqrt(diagj)
                  diag(j) = diagj
                  diagj=1.d0/diagj
               else
                  diag(j)=0.0d0
                  diagj=0.0d0
               endif
            else
               diag(j)=0.0d0
               diagj=0.0d0
            endif
cmae
            ISTRT = XLNZ(J)                                              11
            ISTOP = XLNZ(J+1) - 1                                        11
            IF ( ISTOP .GE. ISTRT )  THEN
               FIRST(J) = ISTRT                                          11
               I = XNZSUB(J)                                             11
               ISUB = NZSUB(I)                                           11
               LINK(J) = LINK(ISUB)                                      12
               LINK(ISUB) = J                                            12
cdir$ ivdep
               DO II = ISTRT, ISTOP
                  ISUB = NZSUB(I)                                        12
cmae              LNZ(II) = ( LNZ(II)-TEMP(ISUB) ) / DIAGJ               12
                  LNZ(II) = ( LNZ(II)-TEMP(ISUB) ) * DIAGJ               12
                  TEMP(ISUB) = 0.0E0                                     12
                  I = I + 1                                              12
               ENDDO
               COUNT = ISTOP - ISTRT + 1                                 12
               OPS  = OPS + COUNT
         ENDIF
      ENDDO
      RETURN
      END                                                                13

c     --------------------
      subroutine figureout (n,avail,iout,st,irank)
c     --------------------
c     mpe
      integer n,avail,irank
      real*8 st(*)
      write(iout,'(a,f25.6)')'               **************'
      write(iout,'(a,f25.6)')'               **** FSPAK ***'
      write(iout,'(a,f25.6)')'               **************'
      write(iout,'(a,f25.6)')'               MPE / IM / MAE'
      write(iout,'(a,f25.6)')'                   Jun 1994'
      write(iout,'(a,f25.6)')
      write(iout,'(a,f25.6)')'              SPARSE STATISTICS'
      write(iout,'(a,I25)')  '      DIMENSION OF MATRIX     =',n
      write(iout,'(a,I25)')  '      RANK                    =',irank
      write(iout,'(a,I25)')  '      STORAGE AVAILABLE       =',avail
      write(iout,'(a,I25)')  '      MAXIMUM NEEDED          ='
     +                                                  ,IVAL(ST(3))
      write(iout,'(a,I25)')  '      NZE IN UPPER TRIANGULAR ='
     +                                                ,IVAL(ST(4))+N
      write(iout,'(a,I25)')  '      NZE IN FACTOR           ='
     +                                                  ,IVAL(ST(2))
      write(iout,'(a,i25)')  '      NO. OF CALLS NUM FACT   ='
     +                                                  ,INT(ST(11))
      write(iout,'(a,i25)')  '      NO. OF CALLS SOLVE      ='
     +                                                  ,INT(ST(12))
      write(iout,'(a,i25)')  '      NO. OF CALLS SPARS SOLV ='
     +                                                  ,INT(ST(13))
      write(iout,'(a,i25)')  '      NO. OF CALLS DET / LDET ='
     +                                                  ,INT(ST(14))
      write(iout,'(a,i25)')  '      NO. OF CALLS SPARS INV  ='
     +                                                  ,INT(ST(15))
      if(st(13).ne.0) then
      write(iout,'(a,i25)')  '      MAX NO. NODES           ='
     +                                                  ,INT(ST(30))
      write(iout,'(a,i25)')  '      MAX PATH LENGTH         ='
     +                                                  ,INT(ST(31))
      write(iout,'(a,f25.3)')'      AVG NO. NODES FPATH     ='
     +                                                ,ST(32)/ST(13)
      write(iout,'(a,f25.3)')'      AVG NO. NODES BPATH     ='
     +                                                ,ST(33)/ST(13)
      write(iout,'(a,f25.3)')'      AVG FPATH LENGTH        ='
     +                                                ,ST(34)/ST(13)
      write(iout,'(a,f25.3)')'      AVG BPATH LENGTH        ='
     +                                                ,ST(35)/ST(13)
      endif
      write(iout,'(a,f25.6)')'      TOTAL CPU TIME IN FSPAK =',st(20)
      write(iout,'(a,f25.6)')'      TIME FOR FINDING ORDER  =',st(21)
      write(iout,'(a,f25.6)')'      TIME FOR SYMBOLIC FAC   =',st(22)
      write(iout,'(a,f25.6)')'      TIME FOR NUMERICAL FAC  =',st(24)
      write(iout,'(a,f25.6)')'      TIME FOR SOLVE          =',st(25)
      write(iout,'(a,f25.6)')'      TIME FOR SPARSE SOLVE   =',st(27)
      write(iout,'(a,f25.6)')'      TIME FOR SPARSE INVERSE =',st(26)
      return
      end

c     -------------
      function ival (st)
c     -------------
      integer st
      ival=st
      return
      end

c     -----------------
      subroutine loadap
c     -----------------
     +                  (n,ia,ja,a,iap,jap,ap,d,ip,tmp)
c this subroutine loads half sparse stored elements of ia,ja,a into half sp
c sparse stored
c structure iap,jap,ap for non-diagonal elements and into d for diagonal el
C elements;
c elements are reordered according to permutation vectors p,ip
c M_P-E Tue Jun 30, 1992 10:23:01
      implicit none
      integer n,ia(1),ja(1),iap(1),jap(1),ip(1),tmp(1),i,j,k,kk
      real*8 a(1),ap(1),d(1)
      do i=1,n
         tmp(i)=0
C        added re: Jan ten Napel suggestion needed for Steve Kachman
C        modifications since checks for zero diagonal
         d(i) = 0
      enddo
c compute no. of entries per upper stored row
      do k=1,n
         i=ip(k)
         do kk=ia(k),ia(k+1)-1
            j=ip(ja(kk))
            if(i.eq.j .and. i.le.n) then
               d(i)=a(kk)
            elseif(i.gt.j .and. i.le.n) then
               tmp(j)=tmp(j)+1
            elseif(i.lt.j .and. j.le.n) then
               tmp(i)=tmp(i)+1
            endif
         enddo
      enddo
c fill iap
      iap(1)=1
      do i=1,n
         iap(i+1)=iap(i)+tmp(i)
         tmp(i)=0
      enddo
c fill jap,ap
      do k=1,n
         i=ip(k)
         do kk=ia(k),ia(k+1)-1
            j=ip(ja(kk))
            if(i.gt.j .and. i.le.n) then
               jap(iap(j)+tmp(j))=i
               ap(iap(j)+tmp(j))=a(kk)
               tmp(j)=tmp(j)+1
            elseif(i.lt.j .and. i.le.n) then
               jap(iap(i)+tmp(i))=j
               ap(iap(i)+tmp(i))=a(kk)
               tmp(i)=tmp(i)+1
            endif
         enddo
      enddo
      return
      end

c     ----------------
      subroutine loadu (n,iap,jap,ap,iu,ju,iju,u,tmp)
c     ----------------
c sparse stored matrix iap,jap,ap is copied into sparse compressed storage
C GE MATRIX IU,JU,IJU,U
c where u allows storage for fill-ins
c M_P-E Tue Jun 30, 1992 11:01:23
      implicit none
      integer n,iap(1),jap(1),iu(1),ju(1),iju(1),i,j,k
      real*8 ap(1),u(1),tmp(1)
      
      do i=1,n
         tmp(i)=0.
      enddo

      do i=1,n
c - -    scatters
         do k=iap(i),iap(i+1)-1
            j=jap(k)
            tmp(j)=ap(k)
         enddo
c - -    loads u and zero out tmp
         do k=iu(i),iu(i+1)-1
            j=ju(iju(i)+k-iu(i))
            u(k)=tmp(j)
            tmp(j)=0.
         enddo
      enddo
      return
      end

c     ----------------
      subroutine lzero (n,l)
c     ----------------
c initialize a logical vector
c MPE
      logical l(1)
      do i=1,n
         l(i)=.false.
      enddo
      return
      end

c     -------------------
      subroutine makepath
c     -------------------
     +           (n,path,iu,ju,iju)
c obtains factorization path stored as a linked list where path(i)
c is next node of the path of element i
c M. Perez-Enciso, Madison,Thu May 21,1992
      integer n,path(n),iu(n+1),ju(1),iju(n)
      do i=1,n
         path(i)=0
         if(iu(i+1).gt.iu(i)) path(i)=ju(iju(i))
      enddo
      return
      end

c     -----------------
      subroutine merlin (n,iu,u,d,tol)
c     -----------------
c MPE
c mae
c skip zero diagonals (machine zero = tol)
      integer n,iu(1)
      real*8 u(1),d(1),tol
      do i=1,n
         if(d(i).gt.tol)then
            do j=iu(i),iu(i+1)-1
               u(j)=-u(j)/d(i)
            enddo
            d(i)=1./(d(i)*d(i))
         endif
      enddo
      return
      end
c     ----------------
      subroutine mkdet (n,d,b,irank,tol,iout,flag)
c     ----------------
c computes determinant, M_P-E
c mae
c skip zero diagonals
      integer i,iout,irank,flag,n,messag
      real*8 b,d(1),tol
      save messag
      data messag /0/
      b=1.
      do i=1,n
         if(d(i).gt.tol)then
            b=b/d(i)
         endif
      enddo
      if(irank.ne.n) then
         flag=54
         if (messag.eq.0) then
            write(iout,*) 'Warning! non full rank matrix'
            write(iout,*) 'Value corresponds to a generalized inverse'
            messag=1
	    print *, 'Sorry, non full rank matrix (--- 1 ---)'
         end if
      endif
      return
      end

c     -----------------
      subroutine mkldet (n,d,b,irank,tol,iout,flag)
c     -----------------
c computes log determinant
c M_P-E Mon Aug 10, 1992 15:43:45
c Mon Feb 28, 1994 10:57:00 (revisited)
c mae
c skip zero diagonals
      integer i,iout,irank,flag,n,messag
      real*8 b,d(1),tol
      save messag
      data messag/0/
      b=0.
      do i=1,n
         if(d(i).gt.tol)then
            b=b-log(d(i))
         endif
      enddo
      if(irank.ne.n) then
         flag=55
         if (messag.eq.0) then
            write(iout,*) 'Warning! non full rank matrix'
            write(iout,*) 'Value corresponds to a generalized inverse'
            messag = 1
            print *, '55: non full rank matrix (--- 2 ---)'
         end if
      endif
      return
      end

c     ------------------
      subroutine permute (n,vec1,vec2,p)
c     ------------------
c permutes vector vec1 according to permutation vector p into vec2
c M_P-E Wed Jul 1, 1992 11:12:01
      integer n,i,p(1)
      real*8 vec1(1),vec2(1)
      do i=1,n
         vec2(p(i))=vec1(i)
      enddo
      return
      end


c     -------------------
      subroutine rowtocol(aia,aja,aa,bia,bja,ba,tmp,n)
c     -------------------
c Matrix a, in sparse i-j-a form and stored row-wise, is converted to b,
c stored column-wise; after the move, the rows of b are sorted in
c increasing order.
c a and b are have n rows and columns, and tmp is
c is a temporary integer vector of size n.
c I. Misztal
      integer n,aia(1),aja(1),bia(1),bja(1),tmp(1),i,j,k
      real*8 aa(1),ba(1)

c count the number of entries in each row of a
      call zero(n,tmp)
      do i=1,n
         do j=aia(i),aia(i+1)-1
            tmp(aja(j))=tmp(aja(j))+1
         enddo
      enddo

c create the row count for b
      bia(1)=1
      do i=1,n
         bia(i+1)=bia(i)+tmp(i)
      enddo

c load a into b
      call zero(n,tmp)
      do i=1,n
         do j=aia(i),aia(i+1)-1
            k=bia(aja(j))+tmp(aja(j))
            bja(k)=i
            ba(k)=aa(j)
            tmp(aja(j))=tmp(aja(j))+1
         enddo
      enddo
      end

c     --------------------
      subroutine rowtocolu (aia,aju,aiju,aa,bia,bja,ba,tmp,n)
c     --------------------
c Matrix a, in sparse i-ju-iju -a form and stored row-wise, is converted to
C TO B,
c stored column-wise; after the move, the rows of b are sorted in
c increasing order.
c a and b are have n rows and columns, and tmp is
c is a temporary integer vector of size n.
c I. Misztal
      integer n,aia(1),bia(1),bja(1),tmp(1),i,j,k,
     *        aju(1),aiju(1)
      real*8 aa(1),ba(1)

c count the number of entries in each row of a
      call zero(n,tmp)
      do i=1,n
         do j=aia(i),aia(i+1)-1
            tmp(aju(aiju(i)-aia(i)+j))=tmp(aju(aiju(i)-aia(i)+j))+1
         enddo
      enddo

c create the row count for b
      bia(1)=1
      do i=1,n
         bia(i+1)=bia(i)+tmp(i)
      enddo

c load a into b
      call zero(n,tmp)
      do i=1,n
         do j=aia(i),aia(i+1)-1
            k=bia(aju(aiju(i)-aia(i)+j))+tmp(aju(aiju(i)-aia(i)+j))
            bja(k)=i
            ba(k)=aa(j)
            tmp(aju(aiju(i)-aia(i)+j))=tmp(aju(aiju(i)-aia(i)+j))+1
         enddo
      enddo
      end

c     ----------------
      subroutine rzero (n,b)
c     ----------------
c initialize a real vector
c MPE
      real*8 b(1)
      do i=1,n
         b(i)=0.
      enddo
      return
      end

c     ------------------
      subroutine segpath
c     ------------------
     +           (n,node,nseg,sstart,spath,path,found)
c computes a segmented path for the nseg nodes stored in node
c M. Perez-Enciso, Madison, Thu May 21, 1992
      implicit none
      integer n,nseg, node(nseg),sstart(nseg+1),spath(1),path(n)
     +        ,ik,k,iseg
      logical found(n)
      ik=0
      do iseg=1,nseg
         k=node(iseg)
         sstart(iseg)=ik+1
         do while(k.ne.0.and..not.found(k))
            ik=ik+1
            spath(ik)=k
            found(k)=.true.
            k=path(k)
         enddo
      enddo
      sstart(nseg+1)=ik+1
c zero out boolean working vector
      do ik=1,sstart(nseg+1)-1
         found(spath(ik))=.false.
      enddo
      return
      end

c     -------------------
      subroutine spermute (n,vect1,vect2,perm)
c     -------------------
c permutes an integer compressed vector
c MPE
      integer vect1(1),vect2(1),perm(1)
      do i=1,n
         vect2(i)=perm(vect1(i))
      enddo
      return
      end

c     -----------------
      subroutine unfold (n,uia,uja,ia,ja,tmp)
c     -----------------
c this subroutine copies a half sparse stored matrix into a sparse
c stored matrix, excluding diagonals
c option=0: no reordering performed, only ia and ja created
c       =1: reordering and sorting performed, a,d also created
c I. Misztal, M. Perez-Enciso
      implicit none
      integer n,uia(1),uja(1),ia(1),ja(1)
     +        ,i,j,k,tmp(1)
c count the number of entries above and below diagonal
      do i=1,n
         tmp(i)=0
      enddo
      do i=1,n
         do k=uia(i),uia(i+1)-1
            j=uja(k)
            if(i.ne.j .and. j.le.n .and. j.ne.0) then
               tmp(i)=tmp(i)+1
               tmp(j)=tmp(j)+1
            endif
         enddo
      enddo
c fill ia ( <--> xadj)
      ia(1)=1
      do i=1,n
         ia(i+1)=ia(i)+tmp(i)
         tmp(i)=0
      enddo
c fill ja ( <--> adjncy )
      do i=1,n
         do k=uia(i),uia(i+1)-1
            j=uja(k)
            if(j.ne.i .and. j.le.n .and. j.ne.0) then
               ja(ia(i)+tmp(i))=j
               ja(ia(j)+tmp(j))=i
               tmp(i)=tmp(i)+1
               tmp(j)=tmp(j)+1
            endif
         enddo
      enddo
      return
      end

c     ---------------
      subroutine zero (n,x)
c     ---------------
      integer n,x(n),i
      do i=1,n
         x(i)=0
      enddo
      end

c     -----------------
      subroutine zerout
c     -----------------
     +                  (b,path,sstart,nseg)
      integer path(1),sstart(1),nseg
      real*8 b(1)
      do i=1,sstart(nseg+1)-1
         b(path(i))=0.
      enddo
      return
      end
C IM 10/24/93-11/23/94
C***********************************************************************
C  Sinverse1 -- inverse OF SPARSE SYMMETRIC POSITIVE DEFINITE SYSTEM OF
C         LINEAR EQUATIONS  Z=inv(M)  GIVEN UT-D-U FACTORIZATION OF M;
C
C         Replaces a in input matrix ia-ja-a with the corresponding
c         inverse elements
C
C         Compared to Sinverse, this version has smaller memory
C         requirements and does not expand ia-ja-a to full storage
C***********************************************************************
      SUBROUTINE  fsinverse1
     *     (N, P, IP, D, IJU,JU,IU,U, TMPu,tmpz,tmpi,tmpi1,
     *      jja,ia,ja,a)
      implicit none
      INTEGER  P(1),  IP(1), IJU(1), JU(1), IU(1), TMPi1(1),
     *           tmpi(1),jja(1),ia(1),ja(1),n,j1,k1
      real*8  D(1), U(1), a(1), TMPu(1),tmpz(1)

      integer i,j,k,iiu
C
C  ADDITIONAL PARAMETERS
C
C    TMPu,tmpz   - real ONE-DIMENSIONAL WORK ARRAYs;  DIMENSION = N
C    jja - integer on-dimensional work array of order ia(n+1)-1
C    tmpi,tmpi1 - integer one-dimensional work array; dimension=N,
c                 tmpi and tmpi1 can occupy the same storage as
c                 tmpu and tmpz
C
C-----------------------------------------------------------------------
      iiu(i,j)=ju(iju(i)-iu(i)+j)
c

c  CAUTION: U is really -U

C----invert in site, z overwrites u
C            tmpu keeps unfolded current column of u;
C            tmpz keeps unfolded current column of z;
      do i=1,n
         tmpz(i)=0
         tmpu(i)=0
      enddo

      do i=n,1,-1
cdir$ ivdep
         do j=iu(i), iu(i+1)-1
            j1=iiu(i,j)
            tmpu(j1)=u(j)
            tmpz(j1)=u(j)*d(j1)
         enddo

c-------off-diagonal elements
         do j=iu(i),iu(i+1)-1
            j1=iiu(i,j)
c           tmpz(j1)=tmpz(j1)+tmpu(j1)*d(j1)
cdir$ ivdep
            do k=iu(j1),iu(j1+1)-1
               k1=iiu(j1,k)
               tmpz(j1)=tmpz(j1)+tmpu(k1)*u(k)
               tmpz(k1)=tmpz(k1)+tmpu(j1)*u(k)
            enddo
          enddo

c------diagonal element last
cdir$ ivdep
         do j=iu(i),iu(i+1)-1
            d(i)=d(i)+u(j)*tmpz(iiu(i,j))
         enddo

c store inverse and zero all nonzeroes in tmp and tmp1
cdir$ ivdep
         do j=iu(i), iu(i+1)-1
            tmpu(iiu(i,j))=0
            u(j)=tmpz(iiu(i,j))
            tmpz(iiu(i,j))=0
         enddo
       enddo

c replace a in ia-ja-a with inverse elements

      do i=1,n
         tmpi1(i)=i
      enddo
c first partly permute ia-ja-a
      call pperm(n,ip,ia,ja,a,tmpi,jja)
      do i=1,n
         do j=iu(ip(i)),iu(ip(i)+1)-1
            tmpu(p(iiu(ip(i),j)))=u(j)
         enddo
         tmpu(i)=d(ip(i))
         do j=ia(i),ia(i+1)-1
            a(j)=tmpu(ja(j))
         enddo
       enddo

c reverse the permutation
      call pperm(n,tmpi1,ia,ja,a,tmpi,jja)
      end


C IM - 10-19-1993-11/23/94
C*********************************************************************
C Partly permute upper-triangular matrix ia-ja-a using permutation vector
C p.
C The permutation is such that it does not change the indices of elements,
C but insures, that after permutation all elements would be in the
C upper diagonal.
C*********************************************************************
        subroutine pperm(n,p,ia,ja,a,tmpia,tmprow)
        implicit none
        integer n,p(1),ia(1),ja(1),tmpia(1),tmprow(1),i,j,k
        real*8 a(1),x

c extra parameters
c tmpia   - integer vector of size n
c tmprow  - integer vector of size ia(n+1)-n

c zero temporary row vector for reorder matrix
        do i=1,n
           tmpia(i)=0
        enddo
c set ja to store column entries and tmprow of row entries
        do i=1,n
           do j=ia(i),ia(i+1)-1
              if (p(i).gt.p(ja(j))) then
                  k=ja(j)
                  ja(j)=i
                  tmprow(j)=k
                 else
                   tmprow(j)=i
                   k=i
                endif
                tmpia(k)=tmpia(k)+1
           enddo
        enddo

c set ia to permuted ia
       do i=1,n
          ia(i+1)=ia(i)+tmpia(i)
          tmpia(i)=ia(i+1)
       enddo

c assign rows to appropriate addresses
        do i=ia(n+1)-1,ia(1),-1
           j=tmprow(i)
           tmpia(j)=tmpia(j)-1
           tmprow(i)=tmpia(j)
        enddo

c final permutation
        do i=1,ia(n+1)-1
10         j=tmprow(i)
           if (i.ne.j) then
c                swap entries i and j
              k=ja(i)
              ja(i)=ja(j)
              ja(j)=k
              k=tmprow(i)
              tmprow(i)=tmprow(j)
              tmprow(j)=k
              x=a(i)
              a(i)=a(j)
              a(j)=x
              goto 10
           endif
        enddo
        end
C  IM - 1/10/92 - 10/24/93
C***********************************************************************
C  Sinverse -- inverse OF SPARSE SYMMETRIC POSITIVE DEFINITE SYSTEM OF
C         LINEAR EQUATIONS  Z=inv(M)  GIVEN UT-D-U FACTORIZATION OF M;
C
C         Expands input matrix ia-ja-a to full storage and repalces
c         cvalues of a with those of the inverse
C***********************************************************************
      SUBROUTINE  fsinverse
     *     (N, P, IP, D, IJU,JU,IU,U, TMP,za,zd,zja,
     *      utia,utja,tmpi,tmpi1,ia,ja,a)
      INTEGER  P(1),  IP(1), IJU(1), JU(1), IU(1), TMPi1(1),
     *          zja(1), tmpi(1),
     *         utia(1),utja(1),ia(1),ja(1)
      real*8  D(1), U(1),  za(1), a(1), TMP(1),zd(1)

      integer i,j,k,l,freeuz,nnzero
C
C  ADDITIONAL PARAMETERS
C
C    TMP   - real ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C    iu,zja,za - inverse, the lower triangular form, in the i-j-a
C                 form; diagonals in ZD
C    utia,utja - U in the i-j-a format, indices only; sorted within rows
C    tmpi,tmpi1 - integer one-dimensional work array; dimension=N
C
C-----------------------------------------------------------------------
C
c  CAUTION: U is really -U
c create indices for U transposed, sorted, za used as temporary
      call rowtocolu(iu,ju,iju,u,utia,utja,za,tmpi,n)

c create I index for Z the same as in U; J index is
c added during computations


C----invert, tmpi contains number of filled elements in each row of z,
C            tmp keeps unfolded current column of z;
c             remember: U is -U
      do i=1,n
         tmpi(i)=0
         tmp(i)=0
      enddo

      do i=n,1,-1
c------scatter existing elements of column i of z into tmp
         do j=iu(i), iu(i)+tmpi(i)-1
            tmp(zja(j))=za(j)
         enddo

c------diagonal element first
         freeuz=iu(i)+tmpi(i)
         tmp(i)=d(i)
         do j=iu(i),iu(i+1)-1
            tmp(i)=tmp(i)+u(j)*tmp(ju(iju(i)-iu(i)+j))
         enddo

         zd(i)=tmp(i)
         nnzero=1
         tmpi1(nnzero)=i

c-------off-diagonal elements
         do j=utia(i+1)-1,utia(i),-1
            freeuz=iu(utja(j))+tmpi(utja(j))
            l=utja(j)
cdir$ ivdep
            do k=iu(l),iu(l+1)-1
               tmp(l)=tmp(l)+u(k)*tmp(ju(iju(l)-iu(l)+k))
            enddo
            za(freeuz)=tmp(l)
            zja(freeuz)=i
            tmpi(utja(j))=tmpi(utja(j))+1
            nnzero=nnzero+1
            tmpi1(nnzero)=l
          enddo

c zero all nonzeroes in tmp
         do j=1,nnzero
            tmp(tmpi1(j))=0
         enddo
       enddo


c in Z, zero all entries not in input matrix A
c ----- since A is upper diagonal, create index structure for A transpose
      call rowtocol(ia,ja,a,utia,utja,a,tmpi,n)

       call zero(n,tmpi)
       call zero(n,tmpi1)
       do i=1,n
c  --------scatter row i of A (upper + lower triangle
          do j=ia(i),ia(i+1)-1
             tmpi(ja(j))=1
          enddo
          do j=utia(i),utia(i+1)-1
             tmpi(utja(j))=1
          enddo
          do j=iu(ip(i)),iu(ip(i)+1)-1
             if (tmpi(p(zja(j))).eq.0) then
                zja(j)=0
             endif
          enddo
c--------zero tmpi and count number of entries in upper+lower diagonal  a
          do j=ia(i),ia(i+1)-1
             tmpi(ja(j))=0
             k=ja(j)
             if (i.ne.k) tmpi1(i)=tmpi1(i)+1
             tmpi1(k)=tmpi1(k)+1
          enddo
          do j=utia(i),utia(i+1)-1
             tmpi(utja(j))=0
          enddo
       enddo


c create row structure for upper+lower diagonal a
      ia(1)=1
      do i=1,n
         ia(i+1)=ia(i)+tmpi1(i)
      enddo

c copy z to a, from half to full stored
      call zero(n,tmpi)
      do i=1,n
         k=p(i)
         a(ia(k)+tmpi(k))=zd(i)
         ja(ia(k)+tmpi(k))=k
         tmpi(k)=tmpi(k)+1
         do j=iu(i),iu(i+1)-1                   
           l=p(i)
           if (zja(j).ne.0) then
                k=p(zja(j)) 
                a(ia(l)+tmpi(l))=za(j)
                a(ia(k)+tmpi(k))=za(j)
                ja(ia(l)+tmpi(l))=k
                ja(ia(k)+tmpi(k))=l
                tmpi(l)=tmpi(l)+1
                tmpi(k)=tmpi(k)+1
            endif
         enddo
      enddo
      return
      end
      subroutine hashmx(y,j1,k1,hashvec1,hashvec2,hashvec3,m,nr)
c stores or acculuates a spare matrix element y(j1,k1) in a hash
c table ind, which rank is m x 3. nr returns the number of stored elemen
c
      integer hashvec1(m),hashvec2(m),j1,k1,m,nr
      real*8 y,hashvec3(m)
      integer*4 iaddr
      integer ie,ieq,izer,iaddress,k,j
      data ie/641/
c
 
      iaddr=433*j1+53*k1
      iaddress=mod(iabs(iaddr),m)+1
 
      do 10 k=1,400
         j=iaddress
         if (hashvec1(iaddress).ne.j1.or.hashvec2(iaddress).ne.k1) then
              ieq=1
            else
              ieq=0
         endif
         if (hashvec1(iaddress).ne.0) then
             izer=1
           else
             izer=0
         endif
 
         if (izer.eq.0 .or. ieq.eq.0) then
             if (izer.eq.0) then
                 hashvec1(iaddress)=j1
                 hashvec2(iaddress)=k1
                 hashvec3(iaddress)=y
                 nr=nr+1
               else
                 hashvec3(iaddress)=hashvec3(iaddress)+y
             endif
c             if (ind(iaddress,3).eq.0) then
c               print*,'hashmx: zero accumulation',j1,k1,ind(iaddress,3)
c             endif
             return
         endif
         iaddress=mod(iaddress+ie-1,m)+1
10    continue
 
      write(*,60)nr,m
60    format(' hash matrix too small,filled',i8,' out of',i8)
      stop
      end

      subroutine hashia(hashvec1,hashvec2,hashvec3,nhash,n,ia,ja,a,m,
     +      tmp)
c copies data form hash to ia-ja-a form, tmp is a temporary array of size n
c Entries are sorted within rows
      integer nhash,n,m,hashvec1(nhash),hashvec2(nhash),ia(n+1),ja(m),
     +      tmp(n)
      real*8 hashvec3(nhash),a(m)
      integer i,imin,imax,j,k,maxcol,maxrow,size

c
c count the number of entries in each column of a
      call zero(n,tmp)
      maxrow=0
      maxcol=0
      do i=1,nhash
         j=hashvec1(i)
         if (j.ne.0) then
            if (j.gt.n) then
                 maxrow=max(maxrow,hashvec1(i))
                 hashvec1(i)=0
            elseif (hashvec2(i).gt.n) then
                 maxcol=max(maxcol,hashvec2(i))
                 hashvec1(i)=0
            else
                 tmp(j)=tmp(j)+1
            endif
         endif
      enddo

c create the row count for b
      ia(1)=1
      do i=1,n
         ia(i+1)=ia(i)+tmp(i)
      enddo

      if (ia(n+1)-1.gt.m) then
         print*,'Too small parameter m in hashia:should be > ', ia(n+1)
         print*,'Increase variable maxnze in param.dat'
         stop
      endif

c load a into b
      call zero(n,tmp)
      do i=1,nhash
         j=hashvec1(i)
         if (j.ne.0) then
            k=ia(j)+tmp(j)
            ja(k)=hashvec2(i)
            a(k)=hashvec3(i)
            tmp(j)=tmp(j)+1
         endif
      enddo
      if ((maxrow.ne.0).or.(maxcol.ne.0)) then
         print*,'HASHIA:declared size=',n,' found columns=',maxcol,
     +         ' and rows =',maxrow
      endif

c sort columns
      do i=1,n
         imin=ia(i)
         imax=ia(i+1)-1
         size=imax-imin+1
         if (size.gt.1) call sortjsp(ja(imin),a(imin),size)
      enddo
      return
      end

      subroutine sortjsp(ja,xa,n)
c sorts integer vector ja in ascending order, and orders the real*8 vector a
c accordingly
      integer i,j,s,l,r,j1,a1,n,x
      integer ja(n),stack(50,2)
      real*8 xa(n),xb

      s=1
      stack(1,1)=1
      stack(1,2)=n
10    l=stack(s,1)
      r=stack(s,2)
      s=s-1
20    i=l
      j=r
      j1=(l+r)/2
      x=ja(j1)

30    if (ja(i).lt.x) then
      i=i+1
      goto 30
      endif

40     if (ja(j).gt.x) then
      j=j-1
      goto 40
      endif
c      print*,i,j,ja(i),ja(j)
      if (i.le.j) then
      a1=ja(i)
      ja(i)=ja(j)
      ja(j)=a1
      xb=xa(i)
      xa(i)=xa(j)
      xa(j)=xb
      i=i+1
      j=j-1
      endif
      if (i.le.j) goto 30
      if (i.lt.r) then
      s=s+1
      stack(s,1)=i
      stack(s,2)=r
      endif
      r=j
      if (l.lt.r) goto 20
      if (s.ne.0) goto 10
      end
C
c      double precision function second()
c      integer sec_100
c      integer first
c      save x,first
c      data first/0/
c
c      if (first.eq.0) then
c        x=sec_100()
c        first=1
c        second=0
c      else
c        second=(sec_100()-x)/100.
c      endif
c      end





C***************************************************************
C***************************************************************
C****     GENMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ****
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE IMPLEMENTS THE MINIMUM DEGREE
C        ALGORITHM.  IT MAKES USE OF THE IMPLICIT REPRESENTATION
C        OF ELIMINATION GRAPHS BY QUOTIENT GRAPHS, AND THE
C        NOTION OF INDISTINGUISHABLE NODES.  IT ALSO IMPLEMENTS
C        THE MODIFICATIONS BY MULTIPLE ELIMINATION AND MINIMUM
C        EXTERNAL DEGREE.
C        ---------------------------------------------
C        CAUTION - THE ADJACENCY VECTOR ADJNCY WILL BE
C        DESTROYED.
C        ---------------------------------------------
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE.
C        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
C        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER
C                 (ANY SMALLER ESTIMATE WILL DO) FOR MARKING
C                 NODES.
C
C     OUTPUT PARAMETERS -
C        PERM   - THE MINIMUM DEGREE ORDERING.
C        INVP   - THE INVERSE OF PERM.
C        NOFSUB - AN UPPER BOUND ON THE NUMBER OF NONZERO
C                 SUBSCRIPTS FOR THE COMPRESSED STORAGE SCHEME.
C
C     WORKING PARAMETERS -
C        DHEAD  - VECTOR FOR HEAD OF DEGREE LISTS.
C        INVP   - USED TEMPORARILY FOR DEGREE FORWARD LINK.
C        PERM   - USED TEMPORARILY FOR DEGREE BACKWARD LINK.
C        QSIZE  - VECTOR FOR SIZE OF SUPERNODES.
C        LLIST  - VECTOR FOR TEMPORARY LINKED LISTS.
C        MARKER - A TEMPORARY MARKER VECTOR.
C
C     PROGRAM SUBROUTINES -
C        MMDELM, MMDINT, MMDNUM, MMDUPD.
C
C***************************************************************
C
      SUBROUTINE  GENMMD ( NEQNS, XADJ, ADJNCY, INVP, PERM,
     1                     DELTA, DHEAD, QSIZE, LLIST, MARKER,
     1                     MAXINT, NOFSUB )
C
C***************************************************************
C
C         INTEGER*2  ADJNCY(1), DHEAD(1) , INVP(1)  , LLIST(1) ,
         INTEGER*4  ADJNCY(1), DHEAD(1) , INVP(1)  , LLIST(1) ,
     1              MARKER(1), PERM(1)  , QSIZE(1)
         INTEGER*4  XADJ(1)
         INTEGER*4  DELTA , EHEAD , I     , MAXINT, MDEG  ,
     1              MDLMT , MDNODE, NEQNS , NEXTMD, NOFSUB,
     1              NUM, TAG
C
C***************************************************************
C
         IF  ( NEQNS .LE. 0 )  RETURN
C
C        ------------------------------------------------
C        INITIALIZATION FOR THE MINIMUM DEGREE ALGORITHM.
C        ------------------------------------------------
         NOFSUB = 0
         CALL  MMDINT ( NEQNS, XADJ, ADJNCY, DHEAD, INVP, PERM,
     1                  QSIZE, LLIST, MARKER )
C
C        ----------------------------------------------
C        NUM COUNTS THE NUMBER OF ORDERED NODES PLUS 1.
C        ----------------------------------------------
         NUM = 1
C
C        -----------------------------
C        ELIMINATE ALL ISOLATED NODES.
C        -----------------------------
         NEXTMD = DHEAD(1)
  100    CONTINUE
             IF  ( NEXTMD .LE. 0 )  GO TO 200
                 MDNODE = NEXTMD
                 NEXTMD = INVP(MDNODE)
                 MARKER(MDNODE) = MAXINT
                 INVP(MDNODE) = - NUM
                 NUM = NUM + 1
                 GO TO 100
C
  200    CONTINUE
C        ----------------------------------------
C        SEARCH FOR NODE OF THE MINIMUM DEGREE.
C        MDEG IS THE CURRENT MINIMUM DEGREE;
C        TAG IS USED TO FACILITATE MARKING NODES.
C        ----------------------------------------
         IF  ( NUM .GT. NEQNS )  GO TO 1000
         TAG = 1
         DHEAD(1) = 0
         MDEG = 2
  300    CONTINUE
             IF  ( DHEAD(MDEG) .GT. 0 )  GO TO 400
                 MDEG = MDEG + 1
                 GO TO 300
  400        CONTINUE
C            -------------------------------------------------
C            USE VALUE OF DELTA TO SET UP MDLMT, WHICH GOVERNS
C            WHEN A DEGREE UPDATE IS TO BE PERFORMED.
C            -------------------------------------------------
             MDLMT = MDEG + DELTA
             EHEAD = 0
C
  500        CONTINUE
                 MDNODE = DHEAD(MDEG)
                 IF  ( MDNODE .GT. 0 )  GO TO 600
                     MDEG = MDEG + 1
                     IF  ( MDEG .GT. MDLMT )  GO TO 900
                         GO TO 500
  600            CONTINUE
C                ----------------------------------------
C                REMOVE MDNODE FROM THE DEGREE STRUCTURE.
C                ----------------------------------------
                 NEXTMD = INVP(MDNODE)
                 DHEAD(MDEG) = NEXTMD
                 IF  ( NEXTMD .GT. 0 )  PERM(NEXTMD) = - MDEG
                 INVP(MDNODE) = - NUM
                 NOFSUB = NOFSUB + MDEG + QSIZE(MDNODE) - 2
                 IF  ( NUM+QSIZE(MDNODE) .GT. NEQNS )  GO TO 1000
C                ----------------------------------------------
C                ELIMINATE MDNODE AND PERFORM QUOTIENT GRAPH
C                TRANSFORMATION.  RESET TAG VALUE IF NECESSARY.
C                ----------------------------------------------
                 TAG = TAG + 1
                 IF  ( TAG .LT. MAXINT )  GO TO 800
                     TAG = 1
                     DO  700  I = 1, NEQNS
                         IF  ( MARKER(I) .LT. MAXINT )  MARKER(I) = 0
  700                CONTINUE
  800            CONTINUE
                 CALL  MMDELM ( MDNODE, XADJ, ADJNCY, DHEAD, INVP,
     1                          PERM, QSIZE, LLIST, MARKER, MAXINT,
     1                          TAG )
                 NUM = NUM + QSIZE(MDNODE)
                 LLIST(MDNODE) = EHEAD
                 EHEAD = MDNODE
                 IF  ( DELTA .GE. 0 )  GO TO 500
  900        CONTINUE
C            -------------------------------------------
C            UPDATE DEGREES OF THE NODES INVOLVED IN THE
C            MINIMUM DEGREE NODES ELIMINATION.
C            -------------------------------------------
             IF  ( NUM .GT. NEQNS )  GO TO 1000
             CALL  MMDUPD ( EHEAD, NEQNS, XADJ, ADJNCY, DELTA, MDEG,
     1                      DHEAD, INVP, PERM, QSIZE, LLIST, MARKER,
     1                      MAXINT, TAG )
             GO TO 300
C
 1000    CONTINUE
         CALL  MMDNUM ( NEQNS, PERM, INVP, QSIZE )
         RETURN
C
      END
C----- SUBROUTINE SMBFCT
C****************************************************************          1.
C****************************************************************          2.
C*********     SMBFCT ..... SYMBOLIC FACTORIZATION       ********          3.
C****************************************************************          4.
C****************************************************************          5.
C                                                                          6.
C     PURPOSE - THIS ROUTINE PERFORMS SYMBOLIC FACTORIZATION               7.
C        ON A PERMUTED LINEAR SYSTEM AND IT ALSO SETS UP THE               8.
C        COMPRESSED DATA STRUCTURE FOR THE SYSTEM.                         9.
C                                                                         10.
C     INPUT PARAMETERS -                                                  11.
C        NEQNS - NUMBER OF EQUATIONS.                                     12.
C        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.                        13.
C        (PERM, INVP) - THE PERMUTATION VECTOR AND ITS INVERSE.           14.
C                                                                         15.
C     UPDATED PARAMETERS -                                                16.
C        MAXSUB - SIZE OF THE SUBSCRIPT ARRAY NZSUB.  ON RETURN,          17.
C               IT CONTAINS THE NUMBER OF SUBSCRIPTS USED                 18.
C                                                                         19.
C     OUTPUT PARAMETERS -                                                 20.
C        XLNZ - INDEX INTO THE NONZERO STORAGE VECTOR LNZ.                21.
C        (XNZSUB, NZSUB) - THE COMPRESSED SUBSCRIPT VECTORS.              22.
C        MAXLNZ - THE NUMBER OF NONZEROS FOUND.                           23.
C        FLAG - ERROR FLAG.  POSITIVE VALUE INDICATES THAT.               24.
C               NZSUB ARRAY IS TOO SMALL.                                 25.
C                                                                         26.
C     WORKING PARAMETERS -                                                27.
C        MRGLNK - A VECTOR OF SIZE NEQNS.  AT THE KTH STEP,               28.
C               MRGLNK(K), MRGLNK(MRGLNK(K)) , .........                  29.
C               IS A LIST CONTAINING ALL THOSE COLUMNS L(*,J)             30.
C               WITH J LESS THAN K, SUCH THAT ITS FIRST OFF-              31.
C               DIAGONAL NONZERO IS L(K,J).  THUS, THE                    32.
C               NONZERO STRUCTURE OF COLUMN L(*,K) CAN BE FOUND           33.
C               BY MERGING THAT OF SUCH COLUMNS L(*,J) WITH               34.
C               THE STRUCTURE OF A(*,K).                                  35.
C        RCHLNK - A VECTOR OF SIZE NEQNS.  IT IS USED TO ACCUMULATE       36.
C               THE STRUCTURE OF EACH COLUMN L(*,K).  AT THE              37.
C               END OF THE KTH STEP,                                      38.
C                   RCHLNK(K), RCHLNK(RCHLNK(K)), ........                39.
C               IS THE LIST OF POSITIONS OF NONZEROS IN COLUMN K          40.
C               OF THE FACTOR L.                                          41.
C        MARKER  - AN INTEGER VECTOR OF LENGTH NEQNS. IT IS USED          42.
C               TO TEST IF MASS SYMBOLIC ELIMINATION CAN BE               43.
C               PERFORMED.  THAT IS, IT IS USED TO CHECK WHETHER          44.
C               THE STRUCTURE OF THE CURRENT COLUMN K BEING               45.
C               PROCESSED IS COMPLETELY DETERMINED BY THE SINGLE          46.
C               COLUMN MRGLNK(K).                                         47.
C                                                                         48.
C****************************************************************         49.
C                                                                         50.
      SUBROUTINE  SMBFCT ( NEQNS, XADJ, ADJNCY, PERM, INVP,               51.
     1                     XLNZ, MAXLNZ, XNZSUB, NZSUB, MAXSUB,           52.
     1                     RCHLNK, MRGLNK, MARKER, FLAG )                 53.
C                                                                         54.
C****************************************************************         55.
C                                                                         56.
         INTEGER ADJNCY(1), INVP(1), MRGLNK(1), NZSUB(1),                 57.
     1           PERM(1), RCHLNK(1), MARKER(1)                            58.
         INTEGER XADJ(1), XLNZ(1), XNZSUB(1),                             59.
     1           FLAG, I, INZ, J, JSTOP, JSTRT, K, KNZ,                   60.
     1           KXSUB, MRGK, LMAX, M, MAXLNZ, MAXSUB,                    61.
     1           NABOR, NEQNS, NODE, NP1, NZBEG, NZEND,                   62.
     1           RCHM, MRKFLG                                             63.
C                                                                         64.
C****************************************************************         65.
C                                                                         66.
C       ------------------                                                67.
C       INITIALIZATION ...                                                68.
C       ------------------                                                69.
        NZBEG = 1                                                         70.
        NZEND = 0                                                         71.
        XLNZ(1) = 1                                                       72.
        DO 100 K = 1, NEQNS                                               73.
           MRGLNK(K) = 0                                                  74.
           MARKER(K) = 0                                                  75.
  100   CONTINUE                                                          76.
C       --------------------------------------------------                77.
C       FOR EACH COLUMN ......... .  KNZ COUNTS THE NUMBER                78.
C       OF NONZEROS IN COLUMN K ACCUMULATED IN RCHLNK.                    79.
C       --------------------------------------------------                80.
        NP1 = NEQNS + 1                                                   81.
        DO 1500  K = 1, NEQNS                                             82.
           KNZ = 0                                                        83.
           MRGK = MRGLNK(K)                                               84.
           MRKFLG = 0                                                     85.
           MARKER(K) = K                                                  86.
           IF (MRGK .NE. 0 ) MARKER(K) = MARKER(MRGK)                     87.
           XNZSUB(K) = NZEND                                              88.
           NODE = PERM(K)                                                 89.
           JSTRT = XADJ(NODE)                                             90.
           JSTOP = XADJ(NODE+1) - 1                                       91.
           IF (JSTRT.GT.JSTOP)  GO TO 1500                                92.
C          -------------------------------------------                    93.
C          USE RCHLNK TO LINK THROUGH THE STRUCTURE OF                    94.
C          A(*,K) BELOW DIAGONAL                                          95.
C          -------------------------------------------                    96.
           RCHLNK(K) = NP1                                                97.
           DO 300 J = JSTRT, JSTOP                                        98.
              NABOR = ADJNCY(J)                                           99.
              NABOR = INVP(NABOR)                                        100.
              IF ( NABOR .LE. K )  GO TO 300                             101.
                 RCHM = K                                                102.
  200            M = RCHM                                                103.
                 RCHM = RCHLNK(M)                                        104.
                 IF ( RCHM .LE. NABOR )  GO TO 200                       105.
                    KNZ = KNZ+1                                          106.
                    RCHLNK(M) = NABOR                                    107.
                    RCHLNK(NABOR) = RCHM                                 108.
                    IF ( MARKER(NABOR) .NE. MARKER(K) )  MRKFLG = 1      109.
  300      CONTINUE                                                      110.
C          --------------------------------------                        111.
C          TEST FOR MASS SYMBOLIC ELIMINATION ...                        112.
C          --------------------------------------                        113.
           LMAX = 0                                                      114.
           IF ( MRKFLG .NE. 0 .OR. MRGK .EQ. 0 ) GO TO 350               115.
           IF ( MRGLNK(MRGK) .NE. 0 ) GO TO 350                          116.
           XNZSUB(K) = XNZSUB(MRGK) + 1                                  117.
           KNZ = XLNZ(MRGK+1) - (XLNZ(MRGK) + 1)                         118.
           GO TO 1400                                                    119.
C          -----------------------------------------------               120.
C          LINK THROUGH EACH COLUMN I THAT AFFECTS L(*,K).               121.
C          -----------------------------------------------               122.
  350      I = K                                                         123.
  400      I = MRGLNK(I)                                                 124.
           IF (I.EQ.0)  GO TO 800                                        125.
              INZ = XLNZ(I+1) - (XLNZ(I)+1)                              126.
              JSTRT = XNZSUB(I) +  1                                     127.
              JSTOP = XNZSUB(I) + INZ                                    128.
              IF (INZ.LE.LMAX)  GO TO 500                                129.
                 LMAX = INZ                                              130.
                 XNZSUB(K) = JSTRT                                       131.
C             -----------------------------------------------            132.
C             MERGE STRUCTURE OF L(*,I) IN NZSUB INTO RCHLNK.            133.
C             -----------------------------------------------            134.
  500         RCHM = K                                                   135.
              DO 700 J = JSTRT, JSTOP                                    136.
                 NABOR = NZSUB(J)                                        137.
  600            M = RCHM                                                138.
                 RCHM = RCHLNK(M)                                        139.
                 IF (RCHM.LT.NABOR)  GO TO 600                           140.
                 IF (RCHM.EQ.NABOR)  GO TO 700                           141.
                    KNZ = KNZ+1                                          142.
                    RCHLNK(M) = NABOR                                    143.
                    RCHLNK(NABOR) = RCHM                                 144.
                    RCHM = NABOR                                         145.
  700         CONTINUE                                                   146.
              GO TO 400                                                  147.
C          ------------------------------------------------------        148.
C          CHECK IF SUBSCRIPTS DUPLICATE THOSE OF ANOTHER COLUMN.        149.
C          ------------------------------------------------------        150.
  800      IF (KNZ.EQ.LMAX)  GO TO 1400                                  151.
C             -----------------------------------------------            152.
C             OR IF TAIL OF K-1ST COLUMN MATCHES HEAD OF KTH.            153.
C             -----------------------------------------------            154.
              IF (NZBEG.GT.NZEND)  GO TO 1200                            155.
                 I = RCHLNK(K)                                           156.
                 DO 900 JSTRT=NZBEG,NZEND                                157.
                    IF (NZSUB(JSTRT)-I)  900, 1000, 1200                 158.
  900            CONTINUE                                                159.
                 GO TO 1200                                              160.
 1000            XNZSUB(K) = JSTRT                                       161.
                 DO 1100 J=JSTRT,NZEND                                   162.
                    IF (NZSUB(J).NE.I)  GO TO 1200                       163.
                    I = RCHLNK(I)                                        164.
                    IF (I.GT.NEQNS)  GO TO 1400                          165.
 1100            CONTINUE                                                166.
                 NZEND = JSTRT - 1                                       167.
C             ----------------------------------------                   168.
C             COPY THE STRUCTURE OF L(*,K) FROM RCHLNK                   169.
C             TO THE DATA STRUCTURE (XNZSUB, NZSUB).                     170.
C             ----------------------------------------                   171.
 1200         NZBEG = NZEND +  1                                         172.
              NZEND = NZEND + KNZ                                        173.
              IF (NZEND.GT.MAXSUB)  GO TO 1600                           174.
              I = K                                                      175.
              DO 1300 J=NZBEG,NZEND                                      176.
                 I = RCHLNK(I)                                           177.
                 NZSUB(J) = I                                            178.
                 MARKER(I) = K                                           179.
 1300         CONTINUE                                                   180.
              XNZSUB(K) = NZBEG                                          181.
              MARKER(K) = K                                              182.
C          --------------------------------------------------------      183.
C          UPDATE THE VECTOR MRGLNK.  NOTE COLUMN L(*,K) JUST FOUND      184.
C          IS REQUIRED TO DETERMINE COLUMN L(*,J), WHERE                 185.
C          L(J,K) IS THE FIRST NONZERO IN L(*,K) BELOW DIAGONAL.         186.
C          --------------------------------------------------------      187.
 1400      IF (KNZ.LE.1)  GO TO 1500                                     188.
              KXSUB = XNZSUB(K)                                          189.
              I = NZSUB(KXSUB)                                           190.
              MRGLNK(K) = MRGLNK(I)                                      191.
              MRGLNK(I) = K                                              192.
 1500      XLNZ(K+1) = XLNZ(K) + KNZ                                     193.
        MAXLNZ = XLNZ(NEQNS) - 1                                         194.
        MAXSUB = XNZSUB(NEQNS)                                           195.
        XNZSUB(NEQNS+1) = XNZSUB(NEQNS)                                  196.
        FLAG = 0                                                         197.
        RETURN                                                           198.
C       ----------------------------------------------------             199.
C       ERROR - INSUFFICIENT STORAGE FOR NONZERO SUBSCRIPTS.             200.
C       ----------------------------------------------------             201.
 1600   FLAG = 1                                                         202.
        RETURN                                                           203.
        END                                                              204.
C***************************************************************
C***************************************************************
C***     MMDINT ..... MULT MINIMUM DEGREE INITIALIZATION     ***
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE PERFORMS INITIALIZATION FOR THE
C        MULTIPLE ELIMINATION VERSION OF THE MINIMUM DEGREE
C        ALGORITHM.
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C        (XADJ,ADJNCY) - ADJACENCY STRUCTURE.
C
C     OUTPUT PARAMETERS -
C        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
C        QSIZE  - SIZE OF SUPERNODE (INITIALIZED TO ONE).
C        LLIST  - LINKED LIST.
C        MARKER - MARKER VECTOR.
C
C***************************************************************
C
      SUBROUTINE  MMDINT ( NEQNS, XADJ, ADJNCY, DHEAD, DFORW,
     1                     DBAKW, QSIZE, LLIST, MARKER )
C
C***************************************************************
C
C         INTEGER*2  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
         INTEGER*4  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
     1              LLIST(1) , MARKER(1), QSIZE(1)
         INTEGER*4  XADJ(1)
         INTEGER*4  FNODE , NDEG  , NEQNS , NODE
C
C***************************************************************
C
         DO  100  NODE = 1, NEQNS
             DHEAD(NODE) = 0
             QSIZE(NODE) = 1
             MARKER(NODE) = 0
             LLIST(NODE) = 0
  100    CONTINUE
C        ------------------------------------------
C        INITIALIZE THE DEGREE DOUBLY LINKED LISTS.
C        ------------------------------------------
         DO  200  NODE = 1, NEQNS
             NDEG = XADJ(NODE+1) - XADJ(NODE) + 1
             FNODE = DHEAD(NDEG)
             DFORW(NODE) = FNODE
             DHEAD(NDEG) = NODE
             IF  ( FNODE .GT. 0 )  DBAKW(FNODE) = NODE
             DBAKW(NODE) = - NDEG
  200    CONTINUE
         RETURN
C
      END
C***************************************************************
C***************************************************************
C**     MMDELM ..... MULTIPLE MINIMUM DEGREE ELIMINATION     ***
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE ELIMINATES THE NODE MDNODE OF
C        MINIMUM DEGREE FROM THE ADJACENCY STRUCTURE, WHICH
C        IS STORED IN THE QUOTIENT GRAPH FORMAT.  IT ALSO
C        TRANSFORMS THE QUOTIENT GRAPH REPRESENTATION OF THE
C        ELIMINATION GRAPH.
C
C     INPUT PARAMETERS -
C        MDNODE - NODE OF MINIMUM DEGREE.
C        MAXINT - ESTIMATE OF MAXIMUM REPRESENTABLE (SHORT)
C                 INTEGER.
C        TAG    - TAG VALUE.
C
C     UPDATED PARAMETERS -
C        (XADJ,ADJNCY) - UPDATED ADJACENCY STRUCTURE.
C        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
C        QSIZE  - SIZE OF SUPERNODE.
C        MARKER - MARKER VECTOR.
C        LLIST  - TEMPORARY LINKED LIST OF ELIMINATED NABORS.
C
C***************************************************************
C
      SUBROUTINE  MMDELM ( MDNODE, XADJ, ADJNCY, DHEAD, DFORW,
     1                     DBAKW, QSIZE, LLIST, MARKER, MAXINT,
     1                     TAG )
C
C***************************************************************
C
C         INTEGER*2  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
         INTEGER*4  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
     1              LLIST(1) , MARKER(1), QSIZE(1)
         INTEGER*4  XADJ(1)
         INTEGER*4  ELMNT , I     , ISTOP , ISTRT , J     ,
     1              JSTOP , JSTRT , LINK  , MAXINT, MDNODE,
     1              NABOR , NODE  , NPV   , NQNBRS, NXNODE,
     1              PVNODE, RLMT  , RLOC  , RNODE , TAG   ,
     1              XQNBR
C
C***************************************************************
C
C        -----------------------------------------------
C        FIND REACHABLE SET AND PLACE IN DATA STRUCTURE.
C        -----------------------------------------------
         MARKER(MDNODE) = TAG
         ISTRT = XADJ(MDNODE)
         ISTOP = XADJ(MDNODE+1) - 1
C        -------------------------------------------------------
C        ELMNT POINTS TO THE BEGINNING OF THE LIST OF ELIMINATED
C        NABORS OF MDNODE, AND RLOC GIVES THE STORAGE LOCATION
C        FOR THE NEXT REACHABLE NODE.
C        -------------------------------------------------------
         ELMNT = 0
         RLOC = ISTRT
         RLMT = ISTOP
         DO  200  I = ISTRT, ISTOP
             NABOR = ADJNCY(I)
             IF  ( NABOR .EQ. 0 )  GO TO 300
                 IF  ( MARKER(NABOR) .GE. TAG )  GO TO 200
                     MARKER(NABOR) = TAG
                     IF  ( DFORW(NABOR) .LT. 0 )  GO TO 100
                         ADJNCY(RLOC) = NABOR
                         RLOC = RLOC + 1
                         GO TO 200
  100                CONTINUE
                     LLIST(NABOR) = ELMNT
                     ELMNT = NABOR
  200    CONTINUE
  300    CONTINUE
C            -----------------------------------------------------
C            MERGE WITH REACHABLE NODES FROM GENERALIZED ELEMENTS.
C            -----------------------------------------------------
             IF  ( ELMNT .LE. 0 )  GO TO 1000
                 ADJNCY(RLMT) = - ELMNT
                 LINK = ELMNT
  400            CONTINUE
                     JSTRT = XADJ(LINK)
                     JSTOP = XADJ(LINK+1) - 1
                     DO  800  J = JSTRT, JSTOP
                         NODE = ADJNCY(J)
                         LINK = - NODE
                         IF  ( NODE )  400, 900, 500
  500                    CONTINUE
                         IF  ( MARKER(NODE) .GE. TAG  .OR.
     1                         DFORW(NODE) .LT. 0 )  GO TO 800
                             MARKER(NODE) = TAG
C                            ---------------------------------
C                            USE STORAGE FROM ELIMINATED NODES
C                            IF NECESSARY.
C                            ---------------------------------
  600                        CONTINUE
                                 IF  ( RLOC .LT. RLMT )  GO TO 700
                                     LINK = - ADJNCY(RLMT)
                                     RLOC = XADJ(LINK)
                                     RLMT = XADJ(LINK+1) - 1
                                     GO TO 600
  700                        CONTINUE
                             ADJNCY(RLOC) = NODE
                             RLOC = RLOC + 1
  800                CONTINUE
  900            CONTINUE
                 ELMNT = LLIST(ELMNT)
                 GO TO 300
 1000    CONTINUE
         IF  ( RLOC .LE. RLMT )  ADJNCY(RLOC) = 0
C        --------------------------------------------------------
C        FOR EACH NODE IN THE REACHABLE SET, DO THE FOLLOWING ...
C        --------------------------------------------------------
         LINK = MDNODE
 1100    CONTINUE
             ISTRT = XADJ(LINK)
             ISTOP = XADJ(LINK+1) - 1
             DO  1700  I = ISTRT, ISTOP
                 RNODE = ADJNCY(I)
                 LINK = - RNODE
                 IF  ( RNODE )  1100, 1800, 1200
 1200            CONTINUE
C                --------------------------------------------
C                IF RNODE IS IN THE DEGREE LIST STRUCTURE ...
C                --------------------------------------------
                 PVNODE = DBAKW(RNODE)
                 IF  ( PVNODE .EQ. 0  .OR.
     1                 PVNODE .EQ. (-MAXINT) )  GO TO 1300
C                    -------------------------------------
C                    THEN REMOVE RNODE FROM THE STRUCTURE.
C                    -------------------------------------
                     NXNODE = DFORW(RNODE)
                     IF  ( NXNODE .GT. 0 )  DBAKW(NXNODE) = PVNODE
                     IF  ( PVNODE .GT. 0 )  DFORW(PVNODE) = NXNODE
                     NPV = - PVNODE
                     IF  ( PVNODE .LT. 0 )  DHEAD(NPV) = NXNODE
 1300            CONTINUE
C                ----------------------------------------
C                PURGE INACTIVE QUOTIENT NABORS OF RNODE.
C                ----------------------------------------
                 JSTRT = XADJ(RNODE)
                 JSTOP = XADJ(RNODE+1) - 1
                 XQNBR = JSTRT
                 DO  1400  J = JSTRT, JSTOP
                     NABOR = ADJNCY(J)
                     IF  ( NABOR .EQ. 0 )  GO TO 1500
                         IF  ( MARKER(NABOR) .GE. TAG )  GO TO 1400
                             ADJNCY(XQNBR) = NABOR
                             XQNBR = XQNBR + 1
 1400            CONTINUE
 1500            CONTINUE
C                ----------------------------------------
C                IF NO ACTIVE NABOR AFTER THE PURGING ...
C                ----------------------------------------
                 NQNBRS = XQNBR - JSTRT
                 IF  ( NQNBRS .GT. 0 )  GO TO 1600
C                    -----------------------------
C                    THEN MERGE RNODE WITH MDNODE.
C                    -----------------------------
                     QSIZE(MDNODE) = QSIZE(MDNODE) + QSIZE(RNODE)
                     QSIZE(RNODE) = 0
                     MARKER(RNODE) = MAXINT
                     DFORW(RNODE) = - MDNODE
                     DBAKW(RNODE) = - MAXINT
                     GO TO 1700
 1600            CONTINUE
C                --------------------------------------
C                ELSE FLAG RNODE FOR DEGREE UPDATE, AND
C                ADD MDNODE AS A NABOR OF RNODE.
C                --------------------------------------
                 DFORW(RNODE) = NQNBRS + 1
                 DBAKW(RNODE) = 0
                 ADJNCY(XQNBR) = MDNODE
                 XQNBR = XQNBR + 1
                 IF  ( XQNBR .LE. JSTOP )  ADJNCY(XQNBR) = 0
C
 1700        CONTINUE
 1800    CONTINUE
         RETURN
C
      END
C***************************************************************
C***************************************************************
C*****     MMDUPD ..... MULTIPLE MINIMUM DEGREE UPDATE     *****
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE UPDATES THE DEGREES OF NODES
C        AFTER A MULTIPLE ELIMINATION STEP.
C
C     INPUT PARAMETERS -
C        EHEAD  - THE BEGINNING OF THE LIST OF ELIMINATED
C                 NODES (I.E., NEWLY FORMED ELEMENTS).
C        NEQNS  - NUMBER OF EQUATIONS.
C        (XADJ,ADJNCY) - ADJACENCY STRUCTURE.
C        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
C        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT)
C                 INTEGER.
C
C     UPDATED PARAMETERS -
C        MDEG   - NEW MINIMUM DEGREE AFTER DEGREE UPDATE.
C        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
C        QSIZE  - SIZE OF SUPERNODE.
C        LLIST  - WORKING LINKED LIST.
C        MARKER - MARKER VECTOR FOR DEGREE UPDATE.
C        TAG    - TAG VALUE.
C
C***************************************************************
C
      SUBROUTINE  MMDUPD ( EHEAD, NEQNS, XADJ, ADJNCY, DELTA,
     1                     MDEG, DHEAD, DFORW, DBAKW, QSIZE,
     1                     LLIST, MARKER, MAXINT, TAG )
C
C***************************************************************
C
C         INTEGER*2  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
         INTEGER*4  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
     1              LLIST(1) , MARKER(1), QSIZE(1)
         INTEGER*4  XADJ(1)
         INTEGER*4  DEG   , DEG0  , DELTA , EHEAD , ELMNT ,
     1              ENODE , FNODE , I     , IQ2   , ISTOP ,
     1              ISTRT , J     , JSTOP , JSTRT , LINK  ,
     1              MAXINT, MDEG  , MDEG0 , MTAG  , NABOR ,
     1              NEQNS , NODE  , Q2HEAD, QXHEAD, TAG
C
C***************************************************************
C
         MDEG0 = MDEG + DELTA
         ELMNT = EHEAD
  100    CONTINUE
C            -------------------------------------------------------
C            FOR EACH OF THE NEWLY FORMED ELEMENT, DO THE FOLLOWING.
C            (RESET TAG VALUE IF NECESSARY.)
C            -------------------------------------------------------
             IF  ( ELMNT .LE. 0 )  RETURN
             MTAG = TAG + MDEG0
             IF  ( MTAG .LT. MAXINT )  GO TO 300
                 TAG = 1
                 DO  200  I = 1, NEQNS
                     IF  ( MARKER(I) .LT. MAXINT )  MARKER(I) = 0
  200            CONTINUE
                 MTAG = TAG + MDEG0
  300        CONTINUE
C            ---------------------------------------------
C            CREATE TWO LINKED LISTS FROM NODES ASSOCIATED
C            WITH ELMNT: ONE WITH TWO NABORS (Q2HEAD) IN
C            ADJACENCY STRUCTURE, AND THE OTHER WITH MORE
C            THAN TWO NABORS (QXHEAD).  ALSO COMPUTE DEG0,
C            NUMBER OF NODES IN THIS ELEMENT.
C            ---------------------------------------------
             Q2HEAD = 0
             QXHEAD = 0
             DEG0 = 0
             LINK = ELMNT
  400        CONTINUE
                 ISTRT = XADJ(LINK)
                 ISTOP = XADJ(LINK+1) - 1
                 DO  700  I = ISTRT, ISTOP
                     ENODE = ADJNCY(I)
                     LINK = - ENODE
                     IF  ( ENODE )  400, 800, 500
C
  500                CONTINUE
                     IF  ( QSIZE(ENODE) .EQ. 0 )  GO TO 700
                         DEG0 = DEG0 + QSIZE(ENODE)
                         MARKER(ENODE) = MTAG
C                        ----------------------------------
C                        IF ENODE REQUIRES A DEGREE UPDATE,
C                        THEN DO THE FOLLOWING.
C                        ----------------------------------
                         IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 700
C                            ---------------------------------------
C                            PLACE EITHER IN QXHEAD OR Q2HEAD LISTS.
C                            ---------------------------------------
                             IF  ( DFORW(ENODE) .EQ. 2 )  GO TO 600
                                 LLIST(ENODE) = QXHEAD
                                 QXHEAD = ENODE
                                 GO TO 700
  600                        CONTINUE
                             LLIST(ENODE) = Q2HEAD
                             Q2HEAD = ENODE
  700            CONTINUE
  800        CONTINUE
C            --------------------------------------------
C            FOR EACH ENODE IN Q2 LIST, DO THE FOLLOWING.
C            --------------------------------------------
             ENODE = Q2HEAD
             IQ2 = 1
  900        CONTINUE
                 IF  ( ENODE .LE. 0 )  GO TO 1500
                 IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 2200
                     TAG = TAG + 1
                     DEG = DEG0
C                    ------------------------------------------
C                    IDENTIFY THE OTHER ADJACENT ELEMENT NABOR.
C                    ------------------------------------------
                     ISTRT = XADJ(ENODE)
                     NABOR = ADJNCY(ISTRT)
                     IF  ( NABOR .EQ. ELMNT )  NABOR = ADJNCY(ISTRT+1)
C                    ------------------------------------------------
C                    IF NABOR IS UNELIMINATED, INCREASE DEGREE COUNT.
C                    ------------------------------------------------
                     LINK = NABOR
                     IF  ( DFORW(NABOR) .LT. 0 )  GO TO 1000
                         DEG = DEG + QSIZE(NABOR)
                         GO TO 2100
 1000                CONTINUE
C                        --------------------------------------------
C                        OTHERWISE, FOR EACH NODE IN THE 2ND ELEMENT,
C                        DO THE FOLLOWING.
C                        --------------------------------------------
                         ISTRT = XADJ(LINK)
                         ISTOP = XADJ(LINK+1) - 1
                         DO  1400  I = ISTRT, ISTOP
                             NODE = ADJNCY(I)
                             LINK = - NODE
                             IF  ( NODE .EQ. ENODE )  GO TO 1400
                             IF  ( NODE )  1000, 2100, 1100
C
 1100                        CONTINUE
                             IF  ( QSIZE(NODE) .EQ. 0 )  GO TO 1400
                             IF  ( MARKER(NODE) .GE. TAG )  GO TO 1200
C                                -------------------------------------
C                                CASE WHEN NODE IS NOT YET CONSIDERED.
C                                -------------------------------------
                                 MARKER(NODE) = TAG
                                 DEG = DEG + QSIZE(NODE)
                                 GO TO 1400
 1200                        CONTINUE
C                            ----------------------------------------
C                            CASE WHEN NODE IS INDISTINGUISHABLE FROM
C                            ENODE.  MERGE THEM INTO A NEW SUPERNODE.
C                            ----------------------------------------
                             IF  ( DBAKW(NODE) .NE. 0 )  GO TO 1400
                             IF  ( DFORW(NODE) .NE. 2 )  GO TO 1300
                                 QSIZE(ENODE) = QSIZE(ENODE) +
     1                                          QSIZE(NODE)
                                 QSIZE(NODE) = 0
                                 MARKER(NODE) = MAXINT
                                 DFORW(NODE) = - ENODE
                                 DBAKW(NODE) = - MAXINT
                                 GO TO 1400
 1300                        CONTINUE
C                            --------------------------------------
C                            CASE WHEN NODE IS OUTMATCHED BY ENODE.
C                            --------------------------------------
                             IF  ( DBAKW(NODE) .EQ.0 )
     1                             DBAKW(NODE) = - MAXINT
 1400                    CONTINUE
                         GO TO 2100
 1500            CONTINUE
C                ------------------------------------------------
C                FOR EACH ENODE IN THE QX LIST, DO THE FOLLOWING.
C                ------------------------------------------------
                 ENODE = QXHEAD
                 IQ2 = 0
 1600            CONTINUE
                     IF  ( ENODE .LE. 0 )  GO TO 2300
                     IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 2200
                         TAG = TAG + 1
                         DEG = DEG0
C                        ---------------------------------
C                        FOR EACH UNMARKED NABOR OF ENODE,
C                        DO THE FOLLOWING.
C                        ---------------------------------
                         ISTRT = XADJ(ENODE)
                         ISTOP = XADJ(ENODE+1) - 1
                         DO  2000  I = ISTRT, ISTOP
                             NABOR = ADJNCY(I)
                             IF  ( NABOR .EQ. 0 )  GO TO 2100
                             IF  ( MARKER(NABOR) .GE. TAG )  GO TO 2000
                                 MARKER(NABOR) = TAG
                                 LINK = NABOR
C                                ------------------------------
C                                IF UNELIMINATED, INCLUDE IT IN
C                                DEG COUNT.
C                                ------------------------------
                                 IF  ( DFORW(NABOR) .LT. 0 )  GO TO 1700
                                     DEG = DEG + QSIZE(NABOR)
                                     GO TO 2000
 1700                            CONTINUE
C                                    -------------------------------
C                                    IF ELIMINATED, INCLUDE UNMARKED
C                                    NODES IN THIS ELEMENT INTO THE
C                                    DEGREE COUNT.
C                                    -------------------------------
                                     JSTRT = XADJ(LINK)
                                     JSTOP = XADJ(LINK+1) - 1
                                     DO  1900  J = JSTRT, JSTOP
                                         NODE = ADJNCY(J)
                                         LINK = - NODE
                                         IF  ( NODE )  1700, 2000, 1800
C
 1800                                    CONTINUE
                                         IF  ( MARKER(NODE) .GE. TAG )
     1                                         GO TO 1900
                                             MARKER(NODE) = TAG
                                             DEG = DEG + QSIZE(NODE)
 1900                                CONTINUE
 2000                    CONTINUE
 2100                CONTINUE
C                    -------------------------------------------
C                    UPDATE EXTERNAL DEGREE OF ENODE IN DEGREE
C                    STRUCTURE, AND MDEG (MIN DEG) IF NECESSARY.
C                    -------------------------------------------
                     DEG = DEG - QSIZE(ENODE) + 1
                     FNODE = DHEAD(DEG)
                     DFORW(ENODE) = FNODE
                     DBAKW(ENODE) = - DEG
                     IF  ( FNODE .GT. 0 )  DBAKW(FNODE) = ENODE
                     DHEAD(DEG) = ENODE
                     IF  ( DEG .LT. MDEG )  MDEG = DEG
 2200                CONTINUE
C                    ----------------------------------
C                    GET NEXT ENODE IN CURRENT ELEMENT.
C                    ----------------------------------
                     ENODE = LLIST(ENODE)
                     IF  ( IQ2 .EQ. 1 )  GO TO 900
                         GO TO 1600
 2300        CONTINUE
C            -----------------------------
C            GET NEXT ELEMENT IN THE LIST.
C            -----------------------------
             TAG = MTAG
             ELMNT = LLIST(ELMNT)
             GO TO 100
C
      END
C***************************************************************
C***************************************************************
C*****     MMDNUM ..... MULTI MINIMUM DEGREE NUMBERING     *****
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE PERFORMS THE FINAL STEP IN
C        PRODUCING THE PERMUTATION AND INVERSE PERMUTATION
C        VECTORS IN THE MULTIPLE ELIMINATION VERSION OF THE
C        MINIMUM DEGREE ORDERING ALGORITHM.
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C        QSIZE  - SIZE OF SUPERNODES AT ELIMINATION.
C
C     UPDATED PARAMETERS -
C        INVP   - INVERSE PERMUTATION VECTOR.  ON INPUT,
C                 IF QSIZE(NODE)=0, THEN NODE HAS BEEN MERGED
C                 INTO THE NODE -INVP(NODE); OTHERWISE,
C                 -INVP(NODE) IS ITS INVERSE LABELLING.
C
C     OUTPUT PARAMETERS -
C        PERM   - THE PERMUTATION VECTOR.
C
C***************************************************************
C
      SUBROUTINE  MMDNUM ( NEQNS, PERM, INVP, QSIZE )
C
C***************************************************************
C
C         INTEGER*2  INVP(1)  , PERM(1)  , QSIZE(1)
         INTEGER*4  INVP(1)  , PERM(1)  , QSIZE(1)
         INTEGER*4  FATHER, NEQNS , NEXTF , NODE  , NQSIZE,
     1              NUM   , ROOT
C
C***************************************************************
C
         DO  100  NODE = 1, NEQNS
             NQSIZE = QSIZE(NODE)
             IF  ( NQSIZE .LE. 0 )  PERM(NODE) = INVP(NODE)
             IF  ( NQSIZE .GT. 0 )  PERM(NODE) = - INVP(NODE)
  100    CONTINUE
C        ------------------------------------------------------
C        FOR EACH NODE WHICH HAS BEEN MERGED, DO THE FOLLOWING.
C        ------------------------------------------------------
         DO  500  NODE = 1, NEQNS
             IF  ( PERM(NODE) .GT. 0 )  GO TO 500
C                -----------------------------------------
C                TRACE THE MERGED TREE UNTIL ONE WHICH HAS
C                NOT BEEN MERGED, CALL IT ROOT.
C                -----------------------------------------
                 FATHER = NODE
  200            CONTINUE
                     IF  ( PERM(FATHER) .GT. 0 )  GO TO 300
                         FATHER = - PERM(FATHER)
                         GO TO 200
  300            CONTINUE
C                -----------------------
C                NUMBER NODE AFTER ROOT.
C                -----------------------
                 ROOT = FATHER
                 NUM = PERM(ROOT) + 1
                 INVP(NODE) = - NUM
                 PERM(ROOT) = NUM
C                ------------------------
C                SHORTEN THE MERGED TREE.
C                ------------------------
                 FATHER = NODE
  400            CONTINUE
                     NEXTF = - PERM(FATHER)
                     IF  ( NEXTF .LE. 0 )  GO TO 500
                         PERM(FATHER) = - ROOT
                         FATHER = NEXTF
                         GO TO 400
  500    CONTINUE
C        ----------------------
C        READY TO COMPUTE PERM.
C        ----------------------
         DO  600  NODE = 1, NEQNS
             NUM = - INVP(NODE)
             INVP(NODE) = NUM
             PERM(NUM) = NODE
  600    CONTINUE
         RETURN
C
      END
