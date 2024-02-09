use cpd_normal;
use std::path::Path;

fn main() {
    println!("Change point detection");
    
    let path_toml = Path::new("test/test_scenario.toml");
    cpd_normal::test_linear(&path_toml, None).unwrap();
}

