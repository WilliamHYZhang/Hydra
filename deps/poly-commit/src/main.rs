use std::env;

use ark_poly::{MultilinearExtension, SparseMultilinearExtension};
use ark_bls12_381::Bls12_381;
use ark_ec::PairingEngine;
use ark_ff::biginteger::BigInteger256;
use ark_std::test_rng;
use ark_poly_commit::multilinear_pc::MultilinearPC;

type E = Bls12_381;
type Fr = <E as PairingEngine>::Fr;

fn convert_bigint(num: u64) -> Fr {
  return Fr::new(BigInteger256([ 064, 0u64, 0u64, num ]));
}

fn main() {
  let args: Vec<String> = env::args().collect();
  println!("{:?}", args);
  let num_gates = args[1].parse::<usize>().unwrap();

  let mut poly_vec = Vec::new();
  for i in 0..num_gates {
    poly_vec.push((i, convert_bigint(args[i+1].parse::<u64>().unwrap())));
  }

  let num_points = args[num_gates+2].parse::<usize>().unwrap();
  let point_dim = args[num_gates+3].parse::<usize>().unwrap();
  let offset = num_gates+4;

  let mut vec_points = Vec::new();
  for i in 0..num_points {
    vec_points.push(vec![]);
    for j in 0..point_dim {
      let index = offset + i * point_dim + j;
      println!("{}", index);
      vec_points[i].push(convert_bigint(args[index].parse::<u64>().unwrap()));
    }
  }

  let poly = SparseMultilinearExtension::from_evaluations(point_dim, &poly_vec);

  let mut rng = test_rng();
  let uni_params = MultilinearPC::<E>::setup(point_dim, &mut rng);

  let nv = poly.num_vars;
  let (ck, vk) = MultilinearPC::<E>::trim(&uni_params, nv);

  let com = MultilinearPC::commit(&ck, &poly);

  for vec_point in vec_points {
    let proof = MultilinearPC::open(&ck, &poly, &vec_point);
    let value = poly.evaluate(&vec_point).unwrap();
    let _result = MultilinearPC::check(&vk, &com, &vec_point, value, &proof);
  }
}
