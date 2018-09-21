!> HLinked List for 2xinteger data TYPE
! + Hash functions

MODULE hash_routines
  USE OpenCMISS
  USE OpenCMISS_Iron
!#include "macros.h"  

!  USE BaseRoutines
!  USE CONSTANTS
!  USE KINDS
!  USE ISO_VARYING_STRING
  implicit none
  
  INTEGER, PARAMETER :: tbl_size = 500 ! a "big" default size for the table

  private  ! by default

  ! TYPEs
  TYPE HLinkedListItem
    integer(CMISSIntg) :: key    
    integer(CMISSIntg) :: value
    TYPE(HLinkedListItem),pointer :: next => NULL()
  END TYPE

  TYPE HLinkedList
    TYPE(HLinkedListItem),pointer :: root => NULL()
    TYPE(HLinkedListItem),pointer :: last => NULL()
  END TYPE

  interface HLinkedList_Add
    module procedure HLinkedList_Add_Data
    module procedure HLinkedList_Add_List
  END interface

! Hash
TYPE HTable
     TYPE(HLinkedList), DIMENSION(:), ALLOCATABLE :: vec
     INTEGER                                      :: vec_len = 0
     LOGICAL                                      :: is_init = .FALSE.
END TYPE HTable

  ! public TYPEs
  public :: HLinkedListItem,HLinkedList, HTable

  ! public subs
  public :: HLinkedList_Add,HLinkedList_Destroy,HLinkedList_Remove_First,HLinkedList_Remove_Last
  public :: HLinkedList_is_Empty,HLinkedList_to_Array

  PUBLIC :: HTable_Init, HTable_Free, HTable_Put, HTable_Get, HTable_Put_All




contains

! -------------------------------------------------------------------

  !> initialises or adds a piece of data to list
  SUBROUTINE HLinkedList_Add_Data(list,key, value,ERR)!,ERROR)

    TYPE(HLinkedList),INTENT(INOUT) :: list
    INTEGER(CMISSIntg),INTENT(IN) :: key, value
    INTEGER(CMISSIntg), INTENT(OUT) :: ERR
!    CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    ! local variables
    TYPE(HLinkedListItem),pointer :: current

    !ENTERS("HLinkedList_Add_Data",ERR,ERROR,*999)

    if (associated(list%root)) then
      current => list%root
       do
         ! key is there already
         IF (current%key==key) then
         EXIT ! do nothing
         END if 

         ! tail is reached
         IF (.NOT. associated(current%next)) THEN
            allocate(current%next)
            current%next%key   = key
            current%next%value = value
            list%last => current%next
         EXIT
         END IF
          
         ! proceed
         current => current%next 

        END DO  
       
     ELSE     ! list is empty: create root 

      allocate(list%root)
      list%root%key   = key
      list%root%value = value
      list%last => list%root
    END IF

    !EXITS("HLinkedList_Add_Data")
    RETURN
!999 ERRORSEXITS("HLinkedList_Add_Data",ERR,ERROR)
!    RETURN 1

  END SUBROUTINE HLinkedList_Add_Data


  SUBROUTINE HLinkedList_Get_Data(list,key, value,ERR)
!,ERROR)

    TYPE(HLinkedList),INTENT(INOUT) :: list
    INTEGER(CMISSIntg),INTENT(IN) :: key
    INTEGER(CMISSIntg),INTENT(OUT) :: value
    INTEGER(CMISSIntg), INTENT(OUT) :: ERR
!    CHARACTER(LEN=*), INTENT(OUT) :: ERROR

    ! local variables
    TYPE(HLinkedListItem),pointer :: current

    ERR = 1000 ! set to 0 if a value is found
    
    if (associated(list%root)) then
      current => list%root
       do
         ! key is found 
         IF (current%key==key) then
         value = current%value 
         ERR = 0

         END if 

         ! tail is reached
         IF (.NOT. associated(current%next)) THEN
         ! value not found!!!
         EXIT
         END IF
          
         ! otherwise proceed
         current => current%next 

        END DO  
    END IF
   ! if list is empty: do nothing 
   
   END SUBROUTINE HLinkedList_Get_Data
! -------------------------------------------------------------------

  !> adds all data from one list (addlist) to another (list)
  Subroutine HLinkedList_Add_List(list,addlist,ERR,ERROR)
    TYPE(HLinkedList),intent(inout) :: list
    TYPE(HLinkedList),intent(in) :: addlist
    INTEGER(CMISSIntg), INTENT(OUT) :: ERR
    CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    ! local variables
    TYPE(HLinkedListItem),pointer :: current

    if (HLinkedList_is_Empty(addlist)) return

    current => addlist%root
    do
      call HLinkedList_Add_Data(list,current%key, current%value, ERR)
     !,ERROR)
      !,*999)
      if (associated(current%next)) then
        current => current%next
      else
        exit
      ENDif
    ENDdo

!    EXITS("HLinkedList_Add_List")
    RETURN
!999 ERRORSEXITS("HLinkedList_Add_List",ERR,ERROR)
!    RETURN 1

  END Subroutine HLinkedList_Add_List

! -------------------------------------------------------------------

  !> removes the first item from list and returns its value in data
  Subroutine HLinkedList_Remove_First(list,key,value,ERR,ERROR,*)
    TYPE(HLinkedList),intent(inout) :: list
    integer(CMISSIntg),intent(out) :: key, value
    INTEGER(CMISSIntg), INTENT(OUT) :: ERR
    CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    ! local variables
    TYPE(HLinkedListItem),pointer :: next
    
    IF (associated(list%root)) THEN
      key   = list%root%key
      value = list%root%value
      next => list%root%next
      deallocate(list%root)
      list%root => next
      IF (associated(list%root)) THEN
        if (.not.associated(list%root%next)) list%last => list%root  ! only one left
      ELSE
        list%last => NULL()
      ENDIF
    ELSE
      write(*,*) ">>> warning: HLinked list is empty and cannot remove first item"
    ENDIF

  END Subroutine HLinkedList_Remove_First

! -------------------------------------------------------------------

  !> removes the last item from list and returns its value in data
  Subroutine HLinkedList_Remove_Last(list,key,value,ERR,ERROR,*)
    TYPE(HLinkedList),intent(inout) :: list
    integer(CMISSIntg),intent(out) :: key,value
    INTEGER(CMISSIntg), INTENT(OUT) :: ERR
    CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    ! local variables
    TYPE(HLinkedListItem),pointer :: current

    if (.not.associated(list%root)) then
      write(*,*) ">>> warning: HLinked list is empty and cannot remove last item"
      return
    ENDif
    current => list%root

    do
      if (associated(current%next)) then
        if (associated(current%next%next)) then
          current => current%next
        else
          ! next one is the last one
          key = current%next%key
          value = current%next%value
          deallocate(current%next)
          current%next => NULL()
          list%last => current
          exit
        ENDif
      else
        ! there must be only one item in the list(?)!
        key   = current%key
        value = current%value
        deallocate(list%root)
        list%root => NULL()
        list%last => NULL()
        exit
      ENDif
    ENDdo

  END Subroutine HLinkedList_Remove_Last

! -------------------------------------------------------------------

  !> will delete and deallocate all items
  Subroutine HLinkedList_Destroy(list,ERR)
!,ERROR)

    TYPE(HLinkedList), INTENT(inout) :: list
    INTEGER(CMISSIntg), INTENT(OUT) :: ERR
 !   CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    ! local variables
    TYPE(HLinkedListItem), POINTER :: current,next

    if (.not.associated(list%root)) return

    current => list%root
    do
      if (associated(current%next)) then
        next => current%next
        deallocate(current)
        current => next
      else
        deallocate(current)
        exit
      ENDif
    ENDdo
    
    list%root => NULL()
    list%last => NULL()

  END SUBROUTINE HLinkedList_Destroy

! -------------------------------------------------------------------

  !> returns true if the list is empty
  Function HLinkedList_is_Empty(list)
    TYPE(HLinkedList),intent(in) :: list
    logical :: HLinkedList_is_Empty

    HLinkedList_is_Empty = .true.
    if (associated(list%root)) HLinkedList_is_Empty = .false.
    
  END Function HLinkedList_is_Empty

! -------------------------------------------------------------------

  !> Returns length of list
  SUBROUTINE HLinkedList_Size(list,n,ERR)
    TYPE(HLinkedList),INTENT(IN) :: list
    INTEGER(CMISSIntg), INTENT(OUT) :: n,ERR
 
    ! local variables
    integer(CMISSIntg) :: i
    TYPE(HLinkedListItem),pointer :: current

    ! return zero-size array if list is empty
    IF (HLinkedList_is_Empty(list)) THEN
      n=0
      RETURN
    END IF

    ! Traversing to find size
    current => list%root
    n=1
    do
      if (associated(current%next)) then
        n=n+1
        current => current%next
      else
        exit
      END IF
    END DO

  END SUBROUTINE HLinkedList_Size

  !> copies out the data to an allocatable array of size 2xn
  Subroutine HLinkedList_to_Array(list,array,ERR,ERROR,*)
    TYPE(HLinkedList),intent(in) :: list
    integer(CMISSIntg),allocatable,intent(out) :: array(:,:)
    INTEGER(CMISSIntg), INTENT(OUT) :: ERR
    CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    ! local variables
    integer(CMISSIntg) :: i,n
    TYPE(HLinkedListItem),pointer :: current

    ! return zero-size array if list is empty
    if (HLinkedList_is_Empty(list)) then
      allocate(array(0,0))
      return
    ENDif

    ! first traversing to find size
    current => list%root
    n=1
    do
      if (associated(current%next)) then
        n=n+1
        current => current%next
      else
        exit
      ENDif
    ENDdo

    ! copy to array
    if (allocated(array)) deallocate(array)
    allocate(array(2,n),stat=err)
    !IF (ERR/=0) CALL ...
    current => list%root
    do i=1,n
      array(1,i)=current%key
      array(2,i)=current%value

      current => current%next
    ENDdo

  END Subroutine HLinkedList_to_Array

  ! initialise hash table
  SUBROUTINE HTable_Init(tbl,ERR, tbl_len)
  !,ERROR)
 
    TYPE (HTable) :: tbl
    INTEGER(CMISSIntg), INTENT(OUT) :: ERR
    !CHARACTER(LEN=*), INTENT(OUT) :: ERROR
    INTEGER(CMISSIntg), OPTIONAL, INTENT(IN)    :: tbl_len

    IF (ALLOCATED(tbl%vec)) DEALLOCATE(tbl%vec)
    IF (PRESENT(tbl_len)) THEN
      ALLOCATE(tbl%vec(0:tbl_len-1))
       tbl%vec_len = tbl_len
    ELSE
       ALLOCATE(tbl%vec(0:tbl_size-1))
       tbl%vec_len = tbl_size
    END IF
   
    tbl%is_init = .TRUE.
  
  END SUBROUTINE HTable_Init

! The first part of the hashing procedure using the string
  ! collating sequence
! NOT NEEDED, not dealing with strings now.
!  ELEMENTAL FUNCTION sum_string(str) RESULT(sig)

  SUBROUTINE HTable_Put(tbl,key,val, ERR)! ,ERROR)
    TYPE (HTable) :: tbl
    INTEGER(CMISSIntg), INTENT(IN)  :: key, val
    INTEGER(CMISSIntg) :: hash
    INTEGER(CMISSIntg), INTENT(OUT) :: ERR
    !CHARACTER(LEN=*),   INTENT(OUT) :: ERROR

    ! compute the index with hash function
    hash = MOD(key,tbl%vec_len) ! which will be replaced by some meaningful function!!

    ! insert data in the list corresponding to hash index
    CALL HLinkedList_Add_Data(tbl%vec(hash),key,val, ERR)
    !, ERROR)

  END SUBROUTINE HTable_Put

  SUBROUTINE HTable_Get(q, n, p, index_found, found, T_key, ERR) ! queries for q and returns index

    INTEGER(CMISSIntg), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: T_key ! the vector of keys
    INTEGER(CMISSIntg), INTENT(IN)  :: q, n, p                         ! queried value, size, and prime number
    INTEGER(CMISSIntg), INTENT(OUT) :: index_found
    LOGICAL, INTENT(OUT) :: found ! if object has been found

    INTEGER(CMISSIntg), INTENT(OUT) :: ERR
    
    INTEGER(CMISSIntg) :: j,block_ind

    found = .FALSE.
    index_found = 0 ! change??? returns k if used inappropriately!!!
    ! copy procedure as in paper
    j = my_alg_module(my_alg_module(T_key(0)*q,p),n)

    block_ind = T_key(j)
    IF (block_ind /= 0) THEN ! if block is non empty
       index_found = my_alg_module(my_alg_module(T_key(block_ind+1)*q,p),T_key(block_ind)**2)+1+block_ind
       IF (T_key(index_found)==q) found = .TRUE.
    END IF

    ! old implementation DELETE ME
    ! compute the index with hash function
    ! hash = MOD(key,tbl%vec_len) ! which will be replaced by some meaningful function!!
    ! insert data in the list corresponding to hash index
    !CALL HLinkedList_Get_Data(tbl%vec(hash),key,val, ERR)
    !, ERROR)

  END SUBROUTINE HTable_Get


  SUBROUTINE HTable_Free(tbl, ERR)
!, ERROR)

    TYPE (HTable) :: tbl
    INTEGER(CMISSIntg)     :: i, low, high
    INTEGER(CMISSIntg), INTENT(OUT) :: ERR
!    CHARACTER(LEN=*),   INTENT(OUT) :: ERROR

    low  = LBOUND(tbl%vec,DIM=1) ! UBOUND returns a scalar corresponding to the upper bound
                                 ! of the array (vec) along dimension 1 (redundant since 1-dimensional)
    high = UBOUND(tbl%vec,DIM=1) 
    IF (ALLOCATED(tbl%vec)) THEN

       DO i=low,high
          ! CALL tbl%vec(i)%free()
          CALL HLinkedList_Destroy(tbl%vec(i),ERR)
!,ERROR)
       END DO

       DEALLOCATE(tbl%vec)

    END IF

    tbl%is_init = .FALSE.
    IF (.NOT.ALLOCATED(tbl%vec)) tbl%vec_len = 0

  END SUBROUTINE HTable_Free

! build representation of table according to Fredman 1984
  SUBROUTINE HTable_Put_All(S_key, S_val, p, T_key, T_val, ERR)! ,ERROR)

    TYPE (HTable) :: W_tbl ! W as a vector of lists
!    TYPE(LinkedList) ::  ! lis
    INTEGER(CMISSIntg), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: T_key ! the hash table
    INTEGER(CMISSIntg), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: T_val
    INTEGER(CMISSIntg), DIMENSION(:), INTENT(IN) :: S_key ! array of keys
    INTEGER(CMISSIntg), DIMENSION(:), INTENT(IN) :: S_val

    INTEGER(CMISSIntg), INTENT(IN) :: p ! the prime number

    INTEGER(CMISSIntg) :: n,i,k,j, card_w, check_sum, check_sum_sq, count_fill, dim_t, prev_ind, index_w, card_w_store
    INTEGER(CMISSIntg), DIMENSION(:), ALLOCATABLE :: ind_in_T
    TYPE(HLinkedListItem) :: W_list_item

    INTEGER(CMISSIntg), INTENT(OUT) :: ERR
    !CHARACTER(LEN=*),   INTENT(OUT) :: ERROR

    n = size(S_key,1)
    ! initialise w as a table of size n
    CALL HTable_Init(W_tbl,ERR, n)



    ! Compute the W_j
    k=2 ! use 2 as start as in paper (but 1 also works for paper case!!)

    DO

     DO i=1,n

      index_w = my_alg_module(my_alg_module((k*S_key(i)),p),n)
      PRINT *, index_w
      CALL HLinkedList_Add_Data(W_tbl%vec(index_w-1),S_key(i),S_val(i), ERR) ! numbered from 0!

         PRINT *, "Success with ", index_w
     END DO

     check_sum_sq =0
     check_sum   = 0
     count_fill  = 0
     DO i=1,n
      CALL HLinkedList_Size(W_tbl%vec(i-1), card_w, ERR)
      IF (card_w /= 0) THEN
      check_sum_sq = check_sum_sq + card_w**2
      check_sum    = check_sum    + card_w
      count_fill   = count_fill +1
      END IF
     END DO

    IF (check_sum_sq<3*n) THEN ! k is correct
     EXIT
    ELSE
     k = k+1 ! try another k
     IF (k>p) RETURN ! ERROR!!!!!!!!
    END IF

   END DO

   ! Allocate T: T(0)=k, T(1),...,T(n) access indices, T(*) |W_j| T(**) k' T(***) = W_j
   dim_T = 1+n+2*count_fill+check_sum_sq
   ALLOCATE (T_key(0:(dim_T-1)),STAT=ERR)
   ALLOCATE (T_val(0:(dim_T-1)),STAT=ERR)
   T_key = 0
   T_val = 0
   T_key(0) = k
   IF (dim_T<=n+1) RETURN ! nothing in here!! 

   prev_ind = n

   DO i=1,n

     ! Cardinality of block W_i 
     CALL HLinkedList_Size(W_tbl%vec(i-1), card_w_store, ERR)

     IF (card_w_store == 0) THEN
       T_key(i) = 0 ! Block is empty
     ELSE

       IF (prev_ind == n) THEN
        T_key(i) = n+1 ! First non-zero entry is just index that follows n
       ELSE 
        ! Compute card of 
        ! CALL HLinkedList_Size(W_tbl%vec(prev_ind-1), card_w, ERR)    
        ! replaced with card_w = card_w_store below!!!!   
        T_key(i) = T_key(prev_ind)+2+card_w**2    
       END IF

       T_key(T_key(i))   = card_w_store ! First  element of each block: card of W
       T_key(T_key(i)+1) = 1            ! Second element of each block: k' (set to 1)

       ALLOCATE (ind_in_T(card_w_store), STAT=ERR)  ! ind_in_T contains the indices where keys will be stored in T for each block    
       W_list_item = W_tbl%vec(i-1)%root            ! take the first item in the block
       j = 1 ! init j for do loop

       DO 
         ! function given by corollary 2
         ! gives LOCAL index
         ! +2 translation
         ind_in_T(j) =  my_alg_module(my_alg_module(T_key(T_key(i)+1)*W_list_item%key,p),card_w_store**2)+2
         IF ((j==1).OR. any(ind_in_T(1:(j-1))/=ind_in_T(j))) THEN ! insert items until conflict

            T_key(ind_in_T(j)+T_key(i)-1) = W_list_item%key 
            T_val(ind_in_T(j)+T_key(i)-1) = W_list_item%value

            IF (associated(W_list_item%next)) W_list_item = W_list_item%next
            
            IF (j==card_w_store) EXIT ! everything filled correctly
            j = j+1  ! next key
            
         ELSE ! conflict: increase k' and start over
             T_key(T_key(i)+1) =T_key(T_key(i)+1)+1 
             T_key(ind_in_T + T_key(i)-1) = 0   ! reset to zero what already filled (keys)
             T_val(ind_in_T + T_key(i)-1) = 0   ! reset to zero what already filled (values)
             ind_in_T = 0                       ! reset to zero indices
             j = 1                              ! reset j
             W_list_item = W_tbl%vec(i-1)%root  ! go back to list start
         END IF
 
       END DO
       
       PRINT *, "IND IN T",  ind_in_T
       DEALLOCATE(ind_in_T)

       prev_ind = i
       card_w = card_w_store
     END IF

   END DO


   print *, "T vector n    ", T_key (1:n)
   print *, "T vector rest ", T_key (n+1:)
   print *, "T value vector n    ", T_val (1:n)
   print *, "T value vector rest ", T_val (n+1:)
   print *, "of size ", size(T_key)

   ! insert data in the list corresponding to hash index
   ! CALL HLinkedList_Add_Data(tbl%vec(hash),key,val, ERR)

  END SUBROUTINE HTable_Put_All

FUNCTION my_alg_module(a,b) RESULT(c)

    INTEGER(CMISSIntg), INTENT(IN)  :: a,b
    INTEGER(CMISSIntg)  :: c

    c = MOD(a,b)
    IF (c==0) c = b

END FUNCTION my_alg_module


END MODULE hash_routines

