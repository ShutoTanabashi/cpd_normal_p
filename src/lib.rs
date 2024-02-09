//  正規分布に従うプロセスに対する変化点検出(Change point detection)・状態追跡法(Process state tracking method)．

pub mod linear;
pub mod step_multiple;

use std::error::Error;
use std::path::Path;

extern crate rayon;
use rayon::prelude::*;

extern crate process_param;
use process_param::{Tau, NumChg};
use process_param::norm::{Parameter, Mu, Sigma2};
extern crate rand_scenario;
use rand_scenario::norm::{RandomScenario, Seed};
extern crate cpd_tools;

/// `cpd_normal`に関するError
#[derive(Debug, Clone)]
pub struct CpdNormError {
    pub message: String,
}

impl std::fmt::Display for CpdNormError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for CpdNormError {
    fn description(&self) -> &str {
        &self.message
    }
}


/// 各期毎の最尤推定量
///
/// 各期毎のサンプルサイズ$ n $が同じデータを計算する．
#[derive(Debug, Clone, PartialEq)]
pub struct MleNorm {
    data: Vec<Vec<f64>>,
    mle: Vec<Parameter>,
    /// 期数$ T $
    pub t: Tau,
    /// 1期毎のサンプルサイズ$ n $ 
    pub n: Tau,
}

impl MleNorm{
    /// データからインスタンスを生成
    ///
    /// # 引数
    /// * `data_norm` -
    /// 最尤推定量を計算する確率変数群．1行毎が各期における取得値（サンプルサイズ$ n $）である．
    ///
    /// # 使用例
    /// ```
    /// # use cpd_normal::MleNorm;
    /// let data = [
    ///     vec!(1.0, 0.0, -1.0),
    ///     vec!(-1.0, 1.0, 0.0),
    ///     vec!(1.0, 0.0, 2.0),
    ///     vec!(2.0, 1.0, 0.0)
    ///     ];
    /// let mle = MleNorm::new(&data).unwrap();
    /// assert_eq!(mle.n, 3);
    /// assert_eq!(mle.t, 4);
    /// ```
    ///
    /// サンプルサイズが期毎に異なると`Err`を返す
    /// ```should_panic
    /// # use cpd_normal::MleNorm;
    /// let data = [
    ///     vec!(1.0, 0.0, -1.0),
    ///     vec!(-1.0, 1.0, 0.0),
    ///     vec!(1.0, 0.0),
    ///     vec!(2.0)
    ///     ];
    /// let mle = MleNorm::new(&data).unwrap(); // panic!
    /// ```
    pub fn new(data_norm: &[Vec<f64>]) -> Result<Self, CpdNormError> {
        let data = data_norm.to_vec();
        let t = data.len() as Tau;
        let n = data[0].len() as Tau;

        // サンプルサイズ確認
        if data[1..].par_iter()
                    .any(|d| (d.len() as Tau) != n) {
            return Err( CpdNormError {
                message: "Sample size is not constant.".to_owned()
            })
        };

        let res_mle = data.par_iter()
                          .map(|d| <Parameter as process_param::Mle>::mle(d) )
                          .collect::<Result<Vec<Parameter>, process_param::ParamError>>();
        let mle;
        match res_mle {
            Err(_) => return Err( CpdNormError {
                message: "The MLE cannot be calculated from the input data.".to_owned()
                      }),
            Ok(m) => mle = m,
        };

        Ok( MleNorm{ data, mle, t, n } )
    }
    
    
    /// 第t期目のデータに対する最尤推定量を取得
    ///
    /// # 引数
    /// * `t` - 期数
    ///
    /// # 使用例
    /// ```
    /// # use cpd_normal::MleNorm;
    /// extern crate process_param;
    /// use process_param::norm::Parameter;
    /// let data = [
    ///     vec!(1.0, 0.0, -1.0),
    ///     vec!(-1.0, 1.0, 0.0),
    ///     vec!(1.0, 0.0, 2.0),
    ///     vec!(3.0, 0.0, 3.0)
    ///     ];
    /// let mle = MleNorm::new(&data).unwrap();
    /// let param_4 = mle.get_mle(&4).unwrap();
    /// assert_eq!(param_4, Parameter::new(2.0, 2.0).unwrap());
    /// ```
    ///
    /// 引数が不適切だと`Err`となる
    /// ```should_panic
    /// # use cpd_normal::MleNorm;
    /// let data = [
    ///     vec!(1.0, 0.0, -1.0),
    ///     vec!(-1.0, 1.0, 0.0),
    ///     vec!(1.0, 0.0, 2.0),
    ///     vec!(3.0, 0.0, 3.0)
    ///     ];
    /// let mle = MleNorm::new(&data).unwrap();
    /// let param_6 = mle.get_mle(&6).unwrap(); // panic!
    /// ```
    pub fn get_mle(&self, t: &Tau) -> Result<Parameter, CpdNormError> {
        match self.check_t(t) {
            Err(e) => Err(e),
            Ok(_) => Ok(self.mle[(*t - 1) as usize].clone()),
        }
    }


    // 期数$ t $が適切な値域かを確認
    //
    // # 引数
    // * `t` - 期数
    fn check_t(&self, t: &Tau) -> Result<(), CpdNormError> {
        // if *t > self.t || *t <= 0 { // DEBUG
        if *t > self.t {
            Err(CpdNormError {
                message: format!("Time step t = {} is out of range: max range is {}", t, self.t)
            })
        } else {
            Ok(())
        }
    }


    // 2個の連続する変化点に対して，値が適切かを確認
    //
    // 値が値域を超える，あるいは前の変化点`t_k_1`の値が後ろの変化点`t_k`の値以上の場合に`Err`となる．
    // * `t_k_1` - 取得する最初の期の一つ手前． `t_k_1` + 1期目からデータを得る．
    // * `t_k` - 取得する最後の期
    fn check_range(&self, t_k_1: &Tau, t_k: &Tau) -> Result<(), CpdNormError> {
        if t_k_1 >= t_k {
            return Err( CpdNormError{
                message: format!(" t_k_1 ({}) needs less than t_k ({})", t_k_1, t_k)
            })
        };
        
        [t_k_1, t_k].iter()
                    .map(|t| self.check_t(t) )
                    .collect::<Result<(), CpdNormError>>()
    }


    /// 第t_k_1 + 1 期から第t_k期までの最尤推定量を取得
    ///
    /// # 引数
    /// * `t_k_1` - 取得する最初の期の一つ手前． `t_k_1` + 1期目からデータを得る．
    /// * `t_k` - 取得する最後の期
    ///
    /// # 使用例（注意点）
    /// 必ず`t_k_1` < `t_k` であること．
    /// 引数の値が不適切だと`Err`を返す．
    /// ```should_panic
    /// # use cpd_normal::MleNorm;
    /// let data = [
    ///     vec!(1.0, 0.0, -1.0),
    ///     vec!(-1.0, 1.0, 0.0),
    ///     vec!(1.0, 0.0, 2.0),
    ///     vec!(3.0, 0.0, 3.0)
    ///     ];
    /// let mle = MleNorm::new(&data).unwrap();
    /// let mle_2_5 = mle.get_mle_range(&2, &5).unwrap(); // panic!
    /// ```
    pub fn get_mle_range(&self, t_k_1: &Tau, t_k: &Tau) -> Result<Vec<Parameter>, CpdNormError> {
        let us_tk1 = *t_k_1 as usize;
        let us_tk = *t_k as usize;
        
        match self.check_range(t_k_1, t_k) {
            Err(e) => Err(e),
            Ok(_) => Ok( self.mle[us_tk1..us_tk].to_vec() )
        }
    }


    /// 第t_k_1 + 1 期から第t_k期までの平均の最尤推定量を取得
    ///
    /// # 引数
    /// * `t_k_1` - 取得する最初の期の一つ手前． `t_k_1` + 1期目からデータを得る．
    /// * `t_k` - 取得する最後の期
    ///
    /// # 使用例（注意点）
    /// 必ず`t_k_1` < `t_k` であること．
    /// 引数の値が不適切だと`Err`を返す．
    pub fn get_mu_range(&self, t_k_1: &Tau, t_k: &Tau) -> Result<Vec<Mu>, CpdNormError> {
        match self.get_mle_range(t_k_1, t_k) {
            Err(e) => Err(e),
            Ok(ps) => Ok(ps.iter()
                           .map(|p| p.mu() )
                           .collect::<Vec<Mu>>()
                        ),
        }
    }


    /// 第t_k_1 + 1 期から第t_k期までの分散の最尤推定量を取得
    ///
    /// # 引数
    /// * `t_k_1` - 取得する最初の期の一つ手前． `t_k_1` + 1期目からデータを得る．
    /// * `t_k` - 取得する最後の期
    ///
    /// # 使用例（注意点）
    /// 必ず`t_k_1` < `t_k` であること．
    /// 引数の値が不適切だと`Err`を返す．
    pub fn get_sigma2_range(&self, t_k_1: &Tau, t_k: &Tau) -> Result<Vec<Sigma2>, CpdNormError> {
        match self.get_mle_range(t_k_1, t_k) {
            Err(e) => Err(e),
            Ok(ps) => Ok(ps.iter()
                           .map(|p| p.sigma2() )
                           .collect::<Vec<Sigma2>>()
                        ),
        }
    }


    /// 最尤推定量をすべて取得
    pub fn mle_all(&self) -> Vec<Parameter> {
        self.mle.clone()
    }
    
    // 第t_k_1 + 1 期から第t_k期までの特性値について，期数の経過で重みづけをした総和 
    // $\sum (t-t_k_1) x_t$
    // を計算
    fn calc_time_weighted_sum(vals: &[f64]) -> f64 {
        vals.iter()
            .enumerate()
            .fold(0.0, |acc, tm| {
                let (i, barx) = tm;
                acc + barx * ((i + 1) as f64)
        })
    }

    /// 第t_k_1 + 1 期から第t_k期までの平均の最尤推定量の総和を取得
    ///
    /// # 引数
    /// * `t_k_1` - 取得する最初の期の一つ手前． `t_k_1` + 1期目からデータを得る．
    /// * `t_k` - 取得する最後の期
    ///
    /// # 使用例（注意点）
    /// 必ず`t_k_1` < `t_k` であること．
    /// 引数の値が不適切だと`Err`を返す．
    pub fn get_mu_sum(&self, t_k_1: &Tau, t_k: &Tau) -> Result<Mu, CpdNormError> {
        let vals = self.get_mu_range(t_k_1, t_k)?;
        let sum = vals.iter()
                      .fold(0.0, |acc, x| acc + x);
        Ok(sum)
    }

    /// 第t_k_1 + 1 期から第t_k期までの平均の最尤推定量について，期数の経過で重みづけをした総和 
    /// $\sum (t-t_k_1) x_t$
    /// を取得
    ///
    /// # 引数
    /// * `t_k_1` - 取得する最初の期の一つ手前． `t_k_1` + 1期目からデータを得る．
    /// * `t_k` - 取得する最後の期
    ///
    /// # 使用例（注意点）
    /// 必ず`t_k_1` < `t_k` であること．
    /// 引数の値が不適切だと`Err`を返す．
    pub fn get_mu_time_weighted_sum(&self, t_k_1: &Tau, t_k: &Tau) -> Result<Mu, CpdNormError> {
        let vals = self.get_mu_range(t_k_1, t_k)?;
        Ok(Self::calc_time_weighted_sum(&vals))
    }

    /// 第 t_k_1 + 1 期から第 t_k 期までの平均の最尤推定量の二乗和
    /// $\sum (x_t^2)$
    /// を取得
    ///
    /// # 引数
    /// * `t_k_1` - 取得する最初の期の一つ手前． `t_k_1` + 1期目からデータを得る．
    /// * `t_k` - 取得する最後の期
    ///
    /// # 使用例（注意点）
    /// 必ず`t_k_1` < `t_k` であること．
    /// 引数の値が不適切だと`Err`を返す．
    pub fn get_mu_square_sum(&self, t_k_1: &Tau, t_k: &Tau) -> Result<Mu, CpdNormError> {
        let vals = self.get_mu_range(t_k_1, t_k)?;
        let ss = vals.iter()
                     .fold(0.0, |acc, v| acc + v.powi(2));
        Ok(ss)
    }

    /// 第t_k_1 + 1 期から第t_k期までの分散の最尤推定量の総和を取得
    ///
    /// # 引数
    /// * `t_k_1` - 取得する最初の期の一つ手前． `t_k_1` + 1期目からデータを得る．
    /// * `t_k` - 取得する最後の期
    ///
    /// # 使用例（注意点）
    /// 必ず`t_k_1` < `t_k` であること．
    /// 引数の値が不適切だと`Err`を返す．
    pub fn get_sigma2_sum(&self, t_k_1: &Tau, t_k: &Tau) -> Result<Sigma2, CpdNormError> {
        let vals = self.get_sigma2_range(t_k_1, t_k)?;
        let sum = vals.iter()
                      .fold(0.0, |acc, x| acc + x);
        Ok(sum)
    }

    /// 第t_k_1 + 1 期から第t_k期までの分散の最尤推定量について，期数の経過で重みづけをした総和 
    /// $\sum (t-t_k_1) x_t$
    /// を取得
    ///
    /// # 引数
    /// * `t_k_1` - 取得する最初の期の一つ手前． `t_k_1` + 1期目からデータを得る．
    /// * `t_k` - 取得する最後の期
    ///
    /// # 使用例（注意点）
    /// 必ず`t_k_1` < `t_k` であること．
    /// 引数の値が不適切だと`Err`を返す．
    pub fn get_sigma2_time_weighted_sum(&self, t_k_1: &Tau, t_k: &Tau) -> Result<Sigma2, CpdNormError> {
        let vals = self.get_sigma2_range(t_k_1, t_k)?;
        Ok(Self::calc_time_weighted_sum(&vals))
    }
}


// 以下シミュレーションを行う関数

// 乱数生成のための補助関数

// TOMLファイルからシナリオを読み取り，それに従う乱数を生成
//
// # 引数
// * `path_scenario` - 乱数生成に用いるシナリオ
// * `seed` - シード値． `None` ならば自動でseed値を生成する．
fn gen_rand_from_tomlfile<P: AsRef<Path>>(path_scenario: &P, seed: Option<Seed>) -> Result<RandomScenario, Box<dyn Error>> {
    let scenario = process_param::norm::Scenario::from_toml(path_scenario)?;
    let rand_from_scenario = match seed {
        Some(s) => RandomScenario::from_scenario_seed(&scenario, s)?, 
        None => RandomScenario::from_scenario(&scenario)?,
    };
    Ok(rand_from_scenario)
}

// TOMLファイルからシナリオを読み取り，それに従う乱数を管理図を併用して生成
//
// # 引数
// * `path_scenario` - 乱数生成に用いるシナリオ
// * `seed` - シード値． `None` ならば自動でseed値を生成する．
fn gen_rand_from_tomlfile_cc<P: AsRef<Path>>(path_scenario: &P, seed: Option<Seed>) -> Result<RandomScenario, Box<dyn Error>> {
    let scenario = process_param::norm::Scenario::from_toml(path_scenario)?;
    let rand_from_scenario = match seed {
        Some(s) => RandomScenario::from_scenario_seed_controlchart(&scenario, s)?, 
        None => RandomScenario::from_scenario_controlchart(&scenario)?,
    };
    println!("Control limits of μ are {:?}", scenario.control_limit_xbar());
    Ok(rand_from_scenario)
}

// ファイル保存のための補助関数

// 乱数列をファイル名を自動生成してTOMLファイルで保存
//
// # 引数
// * `path_scenario` - 乱数生成に用いるシナリオ
// * `rands` - 乱数列
// * `annotation` - ファイル名最後部に付ける注釈．不要ならば `None`．
fn save_rands<P: AsRef<Path>>(path_scenario: &P, rands: &RandomScenario, annotation: Option<&str>) -> Result<(), Box<dyn Error>> {
    let an = match annotation {
        Some(s) => "_".to_owned() + s,
        None => "".to_owned(),
    };
    let str_rand_toml = format!("test/{}_rands{}.toml",
                                        path_scenario.as_ref()
                                                     .file_stem()
                                                     .unwrap()
                                                     .to_str()
                                                     .unwrap(),
                                        an);
    let path_rand_toml = Path::new(&str_rand_toml);
    rands.to_toml(&path_rand_toml)?;
    Ok(())
}


/// ステップ変化が発生する変化点推定をシミュレーション
///
/// 管理図の併用なしバージョン
// 解析には[`step_multiple::CPDMultiStep`]を利用
///
/// # 引数
/// * `path_scenario` - 乱数生成に用いるシナリオ
/// * `seed` - シード値． `None` ならば自動でseed値を生成する．
pub fn test_step_multiple<P: AsRef<Path>>(path_scenario: &P, seed: Option<Seed>) -> Result<(), Box<dyn std::error::Error>> {
    let rands = gen_rand_from_tomlfile(path_scenario, seed)?;
    check_step(path_scenario, &rands)
}

/// 管理図を用いるプロセスにて，ステップ変化が発生する変化点推定をシミュレーション
///
/// 管理図の併用ありバージョン
/// 解析には[`step_multiple::CPDMultiStep`]を利用
///
/// # 引数
/// * `path_scenario` - 乱数生成に用いるシナリオ
/// * `seed` - シード値． `None` ならば自動でseed値を生成する．
pub fn test_step_multiple_cc<P: AsRef<Path>>(path_scenario: &P, seed: Option<Seed>) -> Result<(), Box<dyn std::error::Error>> {
    let rands = gen_rand_from_tomlfile_cc(path_scenario, seed)?;
    check_step(path_scenario, &rands)
}

// 既存の乱数列に対してステップ変化の変化点推定をシミュレーション
//
// 解析には[`step_multiple::CPDMultiStep`]を利用
//
// # 引数
// * `path_rand` - 乱数の格納されたtomlファイル
// * `rands` - 乱数列
fn check_step<P: AsRef<Path>>(path_scenario: &P, rands: &RandomScenario) -> Result<(), Box<dyn Error>> {
    let incontrol = rands.get_init_param();
    let vals = rands.rand_vars();
    let analysis = step_multiple::CPDMultiStep::new(vals, &incontrol)?;

    save_rands(path_scenario, &rands, None)?;

    println!("Model is k={}:\n{:#?}\n", analysis.model().num_change(), analysis.model());

    let path_xlsx = format!("test/{}_analyzed_multistep.xlsx", 
                            path_scenario.as_ref()
                                         .file_stem()
                                         .unwrap()
                                         .to_str()
                                         .unwrap());
    analysis.to_excel(&path_xlsx)?;

    Ok(())
}


/// 線形変化が発生する変化点推定をシミュレーション
///
/// 管理図の併用なしバージョン
/// 解析には[`linear::CPDLinear`]を利用
///
/// # 引数
/// * `path_scenario` - 乱数生成に用いたシナリオ
/// * `seed` - シード値． `None` ならば自動でseed値を生成する．
pub fn test_linear<P: AsRef<Path>>(path_scenario: &P, seed: Option<Seed>) -> Result<(), Box<dyn std::error::Error>> {
    let rands = gen_rand_from_tomlfile(path_scenario, seed)?;
    check_linear(path_scenario, &rands)
}


/// 管理図を用いた場合にて，線形変化が発生する変化点推定をシミュレーション
///
/// 管理図の併用ありバージョン
/// 解析には[`linear::CPDLinear`]を利用
///
/// # 引数
/// * `path_scenario` - 乱数生成に用いたシナリオ
/// * `seed` - シード値． `None` ならば自動でseed値を生成する．
pub fn test_linear_cc<P: AsRef<Path>>(path_scenario: &P, seed: Option<Seed>) -> Result<(), Box<dyn std::error::Error>> {
    let rands = gen_rand_from_tomlfile_cc(path_scenario, seed)?;
    check_linear(path_scenario, &rands)
}

// 既存の乱数をもとに線形変化が発生する場合を想定したシミュレーション
// 
// 解析には[`linear::CPDLinear`]を利用
//
// # 引数
// * `path_scenario` - 乱数生成に用いたシナリオ
// * `rands` - 乱数列
fn check_linear<P: AsRef<Path>>(path_scenario: &P, rands: &RandomScenario) -> Result<(), Box<dyn Error>> {
    let incontrol = rands.get_init_param();
    let vals = rands.rand_vars();

    let analysis = linear::CPDLinear::new(vals, &incontrol)?;

    save_rands(path_scenario, &rands, None)?;

    println!("Model is k={}:\n{:#?}\n", analysis.model().num_change(), analysis.model());

    let path_xlsx = format!("test/{}_analyzed_linear.xlsx", 
                            path_scenario.as_ref()
                                         .file_stem()
                                         .unwrap()
                                         .to_str()
                                         .unwrap());
    analysis.to_excel(&path_xlsx)?;

    Ok(())
}


/// 任意の変化点数を想定して線形変化が発生する変化点推定をシミュレーション
///
/// 管理図の併用なしバージョン
/// 解析には[`linear::CPDLinear`]を利用
///
/// # 引数
/// * `path_scenario` - 乱数生成に用いたシナリオ
/// * `k` - 変化点数
/// * `seed` - シード値． `None` ならば自動でseed値を生成する．
pub fn test_linear_k<P: AsRef<Path>>(path_scenario: &P, k: &NumChg, seed: Option<Seed>) -> Result<(), Box<dyn std::error::Error>> {
    let rands = gen_rand_from_tomlfile(path_scenario, seed)?;
    check_linear_k(path_scenario, &rands, k)
}

/// 任意の変化点数を想定して，管理図を併用した工程にて線形変化が発生する変化点推定をシミュレーション
///
/// 管理図の併用ありバージョン
/// 解析には[`linear::CPDLinear`]を利用
///
/// # 引数
/// * `path_scenario` - 乱数生成に用いたシナリオ
/// * `k` - 変化点数
/// * `seed` - シード値． `None` ならば自動でseed値を生成する．
pub fn test_linear_k_cc<P: AsRef<Path>>(path_scenario: &P, k: &NumChg, seed: Option<Seed>) -> Result<(), Box<dyn std::error::Error>> {
    let rands = gen_rand_from_tomlfile_cc(path_scenario, seed)?;
    check_linear_k(path_scenario, &rands, &k)
}

// 既存の乱数列に対して任意の変化点数を想定した線形変化の変化点推定を実行
//
// 推定には[`linear::CPDLinear`]を利用
//
// # 引数
// * `path_scenario` - 乱数生成に用いたシナリオ
// * `rands` - 乱数列
// * `k` - 変化点数
fn check_linear_k<P: AsRef<Path>>(path_scenario: &P, rands: &RandomScenario, k: &NumChg) -> Result<(), Box<dyn Error>> {
    let incontrol = rands.get_init_param();
    let vals = rands.rand_vars();
    let analysis = linear::CPDLinear::new(vals, &incontrol)?;

    save_rands(path_scenario, &rands, None)?;

    println!("Model is k={}:\n{:#?}\n", analysis.model().num_change(), analysis.model());

    // 最適解の保存
    let path_xlsx_opt = format!("test/{}_analyzed_linear.xlsx", 
                            path_scenario.as_ref()
                                         .file_stem()
                                         .unwrap()
                                         .to_str()
                                         .unwrap());
    analysis.to_excel(&path_xlsx_opt)?;
   
    
    // 任意の変化点個数のモデルを保存
    let path_xlsx_k = format!("test/{}_analyzed_linear_k{}.xlsx", 
                            path_scenario.as_ref()
                                         .file_stem()
                                         .unwrap()
                                         .to_str()
                                         .unwrap(),
                            k);
    analysis.to_excel_k(&path_xlsx_k, k)?;
    
    Ok(())
}


/// 既存の乱数列から，線形変化モデルにおいて尤度を最大にするモデルを出力
///
/// # 引数
/// * `path_rand` - 乱数の格納されたtomlファイル
pub fn test_linear_max_ll<P: AsRef<Path>>(path_rand: &P) -> Result<(), Box<dyn std::error::Error>> {
    let rands = rand_scenario::norm::RandomScenario::from_toml(path_rand)?;
    let analysis = linear::CPDLinear::new(rands.rand_vars(), &rands.get_init_param())?;

    // 最適解の保存
    let path_xlsx_max_ll = format!("test/{}_analyzed_linear_max_ll.xlsx", 
                            path_rand.as_ref()
                                     .file_stem()
                                     .unwrap()
                                     .to_str()
                                     .unwrap());
    analysis.to_excel_max_ll(&path_xlsx_max_ll)?;

    println!("The model given maximum likelihood is saved.");
    Ok(())
}
