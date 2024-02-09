//! 変化点にて線形変化が発生する場合

use super::{CpdNormError, MleNorm};

extern crate simple_excel_writer as xlsx;
use simple_excel_writer::{Row, sheet::ToCellValue};


extern crate process_param;
use process_param::{Tau, NumChg, Process, ChangeType};
use process_param::norm::{Mu, Parameter};
extern crate cpd_tools;
use cpd_tools::dp_tools::calc_dp::{CalcTT, DictTT, CalcDP, order_change_point};
use cpd_tools::dp_tools::CalcDpError;


/// ステップ変化による変化後の平均$ \mu_k $の計算
#[derive(Debug, Clone)]
pub struct Muk {
    muk: Vec<Vec<Mu>>
}


impl Muk {
    /// [`Muk`]インスタンスを作成
    ///
    /// # 引数
    /// * `mle` - 各期毎の最尤推定量
    pub fn new(mle: &MleNorm) -> Result<Self, CpdNormError> {
        let res_muk = <Self as DictTT<Mu, MleNorm>>::calc_value_all(mle, &mle.t);
        
        match res_muk {
            Ok(muk) => Ok( Muk{ muk } ),
            Err(e) => Err( CpdNormError{
                message: format!("Failed to create muk instance because: {e}")
            }),
        }
    }
}


impl CalcTT<Mu, MleNorm> for Muk {
    /// 2個の変化点間の平均を計算
    ///
    /// # 引数
    /// * `data` - 各期毎の最尤推定量
    /// * `t_k_1` - 前の変化点 $t_{k-1}$
    /// * `t_k` - 後ろの変化点 $t_k$
    fn calc_value(data: &MleNorm, t_k_1: Tau, t_k: Tau) -> Result<Mu, CalcDpError> {
        // 最尤推定量の取得
        let mle_mu;
        match data.get_mu_range(&t_k_1, &t_k) {
            Err(e) => return Err( CalcDpError {
                        message: e.message
                      }),
            Ok(ms) => mle_mu = ms,
        };

        let f_tk1 = t_k_1 as f64;
        let f_tk = t_k as f64;
        let gap_tk  = f_tk - f_tk1;

        let sum = mle_mu.iter()
                        .fold(0.0, |acc, x| acc + x);
        Ok( (sum as f64) / gap_tk)
    }
}


impl DictTT<Mu, MleNorm> for Muk {
    fn value_tt_all(&self) -> Vec<Vec<Mu>> {
        self.muk.clone()
    }
}


/// 対数尤度関数（厳密にはモデルに依存する項）
#[derive(Clone)]
pub struct LogLike {
    pub incontorl: Parameter,
    memo: Vec<Vec<Option<(Tau, NumChg, f64)>>>,
    // `Tuple`の中身はそれぞれ`(ひとつ前の変化点，変化点個数， ($ \mu_{k}(\tau_k) $の値， 対数尤度(モデルに依存する値)))`
}


impl LogLike {
    /// [`LogLike`]インスタンスを生成
    ///
    /// # 引数
    /// * `mk` - ステップ変化による平均の推移
    /// * `betak` - 線形変化による平均の傾き
    /// * `mu0` - 管理状態の平均
    pub fn new(mle: &MleNorm, muk: &Muk, incontorl: &Parameter) -> Result<Self, CalcDpError> {
        let memo = Self::calc_memo_all(&(mle.clone(), muk.clone(), incontorl.clone()), &mle.t)?;
        Ok(LogLike{ incontorl: incontorl.clone(), memo })
    }
}


impl CalcTT<f64, (MleNorm, Muk, Parameter)> for LogLike {
    fn calc_value(data: &(MleNorm, Muk, Parameter), t_k_1: Tau, t_k: Tau) -> Result<f64, CalcDpError> {
        // 引数確認
        order_change_point(&t_k_1, &t_k)?;

        let (dict_mle, dict_muk, incontorl) = data;
        let mu0 = incontorl.mu();
        let sigma02 = incontorl.sigma2();

        let muk = dict_muk.value_tt(t_k_1, t_k)?;
        let tgap = (t_k - t_k_1) as f64; 
        let n = dict_mle.n as f64;
       
        let value = n * tgap * (muk - mu0).powi(2) / (2.0 * sigma02);

        Ok(value)
    }
}


impl CalcDP<f64, (MleNorm, Muk, Parameter)> for LogLike {
    fn memo_all(&self) -> Vec<Vec<Option<(Tau, NumChg, f64)>>> {
        self.memo.clone()
    }
}


/// 線形・ステップ変化が混在する状態追跡法
pub struct CPDMultiStep<'t> {
    pub data: &'t [Vec<<Parameter as Process>::Observation>],
    pub incontorl: Parameter,
    #[allow(dead_code)]
    mle: MleNorm,
    #[allow(dead_code)]
    muk: Muk,
    #[allow(dead_code)]
    loglike: LogLike,
    aics: Vec<(Tau, f64)>,
    /// 解析結果
    model: process_param::norm::Scenario,
}


impl<'t> CPDMultiStep<'t> {
    /// データからインスタンスを生成
    ///
    /// # 引数
    /// * `data` - 解析対象のデータ
    /// * `incontorl` - 管理状態のパラメータ
    pub fn new(data: &'t [Vec<<Parameter as Process>::Observation>], incontrol: &Parameter) -> Result<Self, Box<dyn std::error::Error>> {
        let mle = MleNorm::new(data)?;
        let muk = Muk::new(&mle)?;
        let loglike = LogLike::new(&mle, &muk, incontrol)?;
       
        let t = mle.t;
        let (mu0, sigma02) = incontrol.param();

        let aics = Self::calc_aic_all(&loglike, &t)?;

        // AICで最小値を与える変化点個数
        let (k_opt, _) = aics.iter()
                             .reduce(|acc, kv| {
                                if acc.1 <= kv.1 {
                                    acc
                                } else {
                                    kv
                                }
                           }).unwrap();

        let hist_k = loglike.get_value_history(&t, &(*k_opt as Tau))?;
        let mut cps = hist_k.iter()
                               .rev()
                               .map(|h| h.0)
                               .collect::<Vec<Tau>>();
        cps.push(t);

        let mut tk1 = cps[1];
        let mut params = Vec::with_capacity( (*k_opt as usize)+1 );
        let ct_sigma2 = ChangeType::new("Step", &[sigma02])?;

        // 管理状態
        let ct_mu0 = ChangeType::new("Step", &[mu0])?;
        params.push( (tk1, ct_mu0, ct_sigma2) );

        if *k_opt > 0 {
            for tk in cps[2..].iter() {
                let mu_tk = muk.value_tt(tk1, *tk, )?;
                let ct_mu_tk = ChangeType::new("Step", &[mu_tk])?;
                params.push( (*tk, ct_mu_tk, ct_sigma2) );
                tk1 = *tk;
            }
        };

        
        let model = process_param::norm::Scenario::new(mle.n, &params)?;

        Ok( CPDMultiStep{ data, 
                       incontorl: incontrol.clone(),
                       mle,
                       muk, 
                       loglike,
                       aics,
                       model 
        })
    }


    // AICを用いた評価関数
    fn calc_aic(loglike: &LogLike, t: &Tau, k: &NumChg) -> Result<f64, CalcDpError> {
       let ll = loglike.get_value(t, k)?;
       let bk = if *k * 2 < *t {
                    (*k as f64) * 2.0
                } else {
                    *t as f64
                };
       Ok( -2.0 * ll + 2.0 * bk )
    }


    // すべての変化点個数に対してAICを計算
    fn calc_aic_all(loglike: &LogLike, t: &Tau) -> Result<Vec<(Tau, f64)>, CalcDpError> {
        let max_k = *t - 1;
        let aics = (0..=max_k).map(|k| Self::calc_aic(loglike, t, &k) )
                              .collect::<Result<Vec<f64>, CalcDpError>>()?;
        let t_val = aics.iter()
                        .enumerate()
                        .map(|(t,v)| (t as Tau, *v))
                        .collect();
        Ok(t_val)
    }


    /// 推定モデルを取得
    pub fn model(&self) -> process_param::norm::Scenario {
        self.model.clone()
    }


    /// 変化点個数毎のAICを返す
    pub fn k_aic(&self) -> Vec<(Tau, f64)> {
        self.aics.clone()
    }



    /// 解析結果をExcelで出力
    ///
    /// # 引数
    /// * `xlsx_path` - 保存先のExcelファイルパス
    pub fn to_excel(&self, xlsx_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut wb = xlsx::Workbook::create(xlsx_path);
        
        let mut sheet_mle = wb.create_sheet("Mu");
        { 
            let model_params_all = self.model().decomplession()?;
            let mle_all = self.mle.mle_all();

            wb.write_sheet(&mut sheet_mle, |sw| {
                sw.append_row(xlsx::row!["TimeStep(t)", "Average of Data", "Estimated values of mu", (), "Variance of data"])?;
                for (i, mpdt) in model_params_all.iter().zip(mle_all.iter()).enumerate() {
                    let (mp, dt) = mpdt;
                    let hat_mu = mp.mu();
                    let d_ave = dt.mu();
                    let d_var = dt.sigma2();
                    let timestep = (i + 1) as f64;
                    sw.append_row(xlsx::row![timestep.to_cell_value(), d_ave.to_cell_value(), hat_mu.to_cell_value(), (), d_var.to_cell_value()])?;
                }
                Ok(())
            })?;
        }

        let mut sheet_model = wb.create_sheet("Model");
        {
            wb.write_sheet(&mut sheet_model, |sw| {
                sw.append_row(xlsx::row!["k: Position of change points", "Change point","μ_k"])?;
                for (k, mk) in self.model().parameters().iter().enumerate() {
                    let tk = mk.tau() as f64;
                    let m = mk.mu().get_setting()[0];
                    sw.append_row(xlsx::row![(k as f64).to_cell_value(), tk.to_cell_value(), (), (m as f64).to_cell_value()])?;
                }
                Ok(())
            })?;
        }

        let mut sheet_aic = wb.create_sheet("AIC");
        {
            wb.write_sheet(&mut sheet_aic, |sw| {
                sw.append_row(xlsx::row!["K: Number of change points", "Evaluated value (=AIC/2-Const)"])?;
                for (k, aic) in self.k_aic().iter() {
                    sw.append_row(xlsx::row![(*k as f64).to_cell_value(), aic.to_cell_value()])?;
                }
                Ok(())
            })?;
        }


        let mut sheet_data = wb.create_sheet("SampleData");
        {
            wb.write_sheet(&mut sheet_data, |sw| {
                sw.append_row(xlsx::row!["time step"])?;
                for (i, di) in self.data.iter().enumerate() {
                    let mut di_mut = di.clone();
                    let mut d_vec = vec![i as f64];
                    d_vec.append(&mut di_mut);
                    let d_row = simple_excel_writer::sheet::Row::from_iter(d_vec.into_iter());
                    sw.append_row(d_row)?;
                }
                Ok(())
            })?;
        }

        Ok(())
    }
}
