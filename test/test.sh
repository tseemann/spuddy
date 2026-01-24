#!/usr/bin/env bats

setup () {
  name="spuddy"
  bats_require_minimum_version 1.5.0
  dir=$(dirname "$BATS_TEST_FILENAME")
  cd "$dir"
  exe="$dir/../bin/$name"
  fasta="JKD6159.fna.gz"
}

@test "Script syntax check" {
  run -0 perl -c "$exe"
}
@test "Version" {
  run -0 $exe -v
  [[ "$output" =~ "$name " ]]
}
@test "Help" {
  run -0 $exe -h
  [[ "$output" =~ "species" ]]
}
@test "Ciiation" {
  run -0 $exe -C
  [[ "$output" =~ "Seemann" ]]
}
@test "No parameters" {
  run ! $exe
}
@test "Bad option" {
  run ! $exe -Q
}
@test "Secret MOTD dump" {
  run $exe -M
  [[ "$output" =~ "potato" ]]
}
@test "Missing database" {
  run ! $exe -D does_not_exist $fasta
  [[ "$output" =~ "ERROR:" ]]
}
@test "Incorrectly named database" {
  run ! $exe -D $exe $fasta
  [[ "$output" =~ "ERROR:" ]]
}
@test "Illegal -p parameter (float)" {
  run ! $exe -p 3.1415 $fasta
  [[ "$output" =~ "ERROR:" ]]
}
@test "Illegal -c parameter (int)" {
  run ! $exe -c 0 $fasta
  [[ "$output" =~ "ERROR:" ]]
}
