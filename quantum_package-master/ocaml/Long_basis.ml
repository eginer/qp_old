open Core.Std;;
open Qptypes;;

type t = (Symmetry.Xyz.t * Gto.t * Nucl_number.t ) list with sexp

let of_basis b =
  let rec do_work accu = function
    | [] -> accu
    | (g,n)::tail ->
        begin
          let new_accu = 
            Symmetry.Xyz.of_symmetry g.Gto.sym 
            |> List.map ~f:(fun x-> (x,g,n)) 
          in
          do_work (new_accu@accu) tail
        end
  in
  do_work [] b
  |> List.rev
;;

let to_basis b =
  let rec do_work accu = function
  | [] -> List.rev accu
  | (s,g,n)::tail -> 
    let first_sym = 
      Symmetry.Xyz.of_symmetry g.Gto.sym
      |> List.hd_exn
    in
    let new_accu = 
      if ( s = first_sym ) then
        (g,n)::accu
      else
        accu
    in
    do_work new_accu tail
  in
  do_work [] b 
;;

let to_string b =
  let middle = List.map ~f:(fun (x,y,z) -> 
     "( "^((Int.to_string (Nucl_number.to_int z)))^", "^
     (Symmetry.Xyz.to_string x)^", "^(Gto.to_string y)
     ^" )"
  ) b
  |> String.concat ~sep:",\n"
  in "("^middle^")"
;;

include To_md5;;
let to_md5 = to_md5 sexp_of_t
;;

