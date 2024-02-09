//! 変化点にて線形変化が発生する場合

use super::{CpdNormError, MleNorm};

extern crate simple_excel_writer as xlsx;
use simple_excel_writer::{Row, sheet::ToCellValue};


extern crate process_param;
use process_param::{Tau, NumChg, Process, ChangeType};
use process_param::norm::{Mu, Sigma2, Parameter};
extern crate cpd_tools;
use cpd_tools::dp_tools::calc_dp::{CalcTT, DictTT, CalcDPWithVari, order_change_point};
use cpd_tools::dp_tools::CalcDpError;


/// 線形変化による平均の傾き$ \beta_k $の計算
#[derive(Debug, Clone)]
pub struct Betak {
    betak: Vec<Vec<Mu>>
}


impl Betak {
    /// [`Betak`]インスタンスを作成
    ///
    /// # 引数
    /// * `mle` - 各期毎の最尤推定量
    pub fn new(mle: &MleNorm) -> Result<Self, CpdNormError> {
        let res_betak = <Self as DictTT<Mu, MleNorm>>::calc_value_all(mle, &mle.t);
        
        match res_betak {
            Ok(betak) => Ok( Betak{ betak } ),
            Err(e) => Err( CpdNormError{
                message: format!("Failed to create batak instance because: {e}")
            }),
        }
    }


    /// 変化点における線形関数の接合部分$ \mu_{k-1} (\tau_{k-1}) = \mu_k (\tau_{k-1}) $の値をもとにして，平均$ \mu_k $の傾きを計算
    /// 
    /// # 引数
    /// * `t_k_1` - 前の変化点 $t_{k-1}$
    /// * `t_k` - 後ろの変化点 $t_k$
    /// * `mu_tk1` - t_k_1における平均の理論値
    pub fn get_beta_k(&self, t_k_1: &Tau, t_k: &Tau, mu_tk1: &Mu) -> Result<Mu, CalcDpError> {
        order_change_point(t_k_1, t_k)?;

        let m_tk = *t_k as Mu;
        let m_tk1 = *t_k_1 as Mu;
        let gap_tk = m_tk - m_tk1;
        
        let term_cum = self.value_tt(*t_k_1, *t_k)?;
        let term_mk1 = (*mu_tk1 * 3.0) / (2.0 * gap_tk + 1.0);
        Ok( term_cum - term_mk1)
    }
}


impl CalcTT<Mu, MleNorm> for Betak {
    /// 2個の変化点間の線形変化による平均の増分$ beta_k $のうち，接合点$ \mu_{k-1}  $に依存しない項を計算
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
        
        // 平均の最尤推定量$ \bar{X}_t $に関する総和の項
        let cums = mle_mu.iter()
                         .enumerate()
                         .fold(0.0, |acc, tm| {
            let (i, barx) = tm;
            acc + barx * ((i + 1) as f64)
        });

        let coef_denom = gap_tk * (gap_tk + 1.0) * (2.0 * gap_tk + 1.0); 
        Ok( 6.0 * cums / coef_denom )
    }
}


impl DictTT<Mu, MleNorm> for Betak {
    fn value_tt_all(&self) -> Vec<Vec<Mu>> {
        self.betak.clone()
    }
}


/// 対数尤度関数（厳密にはモデルに依存する項）
#[derive(Clone)]
pub struct LogLike {
    pub mu0: Mu,
    memo: Vec<Vec<Option<(Tau, NumChg, Mu, f64)>>>,
    // `Tuple`の中身はそれぞれ`(ひとつ前の変化点，変化点個数， ($ \mu_{k}(\tau_k) $の値， 対数尤度(モデルに依存する値)))`
}


impl LogLike {
    /// [`LogLike`]インスタンスを生成
    ///
    /// # 引数
    /// * `mle` - ステップ変化による平均の推移
    /// * `betak` - 線形変化による平均の傾き
    /// * `mu0` - 管理状態の平均
    pub fn new(mle: &MleNorm, betak: &Betak, mu0: Mu) -> Result<Self, CalcDpError> {
        let memo = Self::calc_memo_all(&(mle.clone(), betak.clone(), mu0), &mle.t)?;
        Ok(LogLike{ mu0, memo })
    }
}


impl CalcDPWithVari<f64, Mu, (MleNorm, Betak, Mu)> for LogLike {
    fn calc_value(data: &(MleNorm, Betak, Mu), t_k_1: &Tau, t_k: &Tau, vari_tk1: &Mu) -> Result<(Mu, f64), CalcDpError> {
        // 引数確認
        order_change_point(&t_k_1, &t_k)?;

        let (dict_mle, dict_betak, mu0) = data;

        let betak = dict_betak.get_beta_k(&t_k_1, &t_k, &vari_tk1)?;
        let tgap = (*t_k - *t_k_1) as f64; 

        let mles;
        match dict_mle.get_mu_range(&t_k_1, &t_k) {
            Ok(m) => mles = m,
            Err(e) => {
                return Err( CalcDpError{
                    message: e.message
                });
            },
        };
        

        let cnst_term = tgap * (
                            vari_tk1.powi(2)
                            +
                            betak * (tgap + 1.0) * (
                                *vari_tk1
                                +
                                betak * (2.0 * tgap + 1.0) / 6.0
                            )
                        );
        let sum = mles.iter()
                      .enumerate()
                      .fold(0.0, |acc, im| {
                          let (i, mle) = im;
                          acc + mle * (
                              *vari_tk1 - mu0 + betak * (1.0 + i as f64)
                          )
                      });
        let ll_tk = 2.0 * sum - cnst_term;

        let mu_tk = *vari_tk1 + betak * tgap;
        
        Ok( (mu_tk, ll_tk) )
    }

    
    fn calc_value_terminal(data: &(MleNorm, Betak, Mu), t_k: &Tau) -> Result<(Mu, f64), CalcDpError> {
        let mu0 = data.2;
        let val = - 1.0 * mu0.powi(2) * (*t_k as f64);
        Ok( (mu0, val ) )
    }


    fn memo_all(&self) -> Vec<Vec<Option<(Tau, NumChg, Mu, f64)>>> {
        self.memo.clone()
    }
}


/// 線形・ステップ変化が混在する状態追跡法
pub struct CPDLinear<'t> {
    pub data: &'t [Vec<<Parameter as Process>::Observation>],
    pub incontrol: Parameter,
    #[allow(dead_code)]
    mle: MleNorm,
    #[allow(dead_code)]
    betak: Betak,
    #[allow(dead_code)]
    loglike: LogLike,
    aics: Vec<(Tau, f64)>,
    /// 解析結果
    model: process_param::norm::Scenario,
}


impl<'t> CPDLinear<'t> {
    /// データからインスタンスを生成
    ///
    /// # 引数
    /// * `data` - 解析対象のデータ
    /// * `incontrol` - 管理状態のパラメータ
    pub fn new(data: &'t [Vec<<Parameter as Process>::Observation>], incontrol: &Parameter) -> Result<Self, Box<dyn std::error::Error>> {
        let mle = MleNorm::new(data)?;
        let betak = Betak::new(&mle)?;
        let loglike = LogLike::new(&mle, &betak, incontrol.mu())?;
       
        let t = mle.t;
        let n = mle.n;
        let (_, sigma02) = incontrol.param();

        let aics = Self::calc_aic_all(&loglike, &n, &t, &sigma02)?;

        // AICで最小値を与える変化点個数
        let (k_opt, _) = aics.iter()
                             .reduce(|acc, kv| {
                                if acc.1 <= kv.1 {
                                    acc
                                } else {
                                    kv
                                }
                           }).unwrap();

        let model = Self::get_model(&mle.t, k_opt, incontrol, &mle, &betak, &loglike)?;

        Ok( CPDLinear{ data, 
                       incontrol: incontrol.clone(),
                       mle,
                       betak, 
                       loglike,
                       aics,
                       model 
        })
    }


    // 任意の変化点個数における最適なモデルを取得
    // 内部関数版
    //
    // # 引数
    // * `t` - 期数
    // * `k` - 変化点個数
    // * `incontrol` - 管理状態のパラメータ
    // * `mle` - サンプル平均
    // * `betak` - 傾きの最尤推定量
    // * `loglike` - 最尤推定量
    fn get_model(t: &Tau, k: &NumChg, incontrol: &Parameter, mle: &MleNorm, betak: &Betak, loglike: &LogLike) -> Result<process_param::norm::Scenario, Box<dyn std::error::Error>> {
        let (mu0, sigma02) = incontrol.param();
        
        let hist_k = loglike.get_value_history(t, k)?;
        let mut cps = hist_k.iter()
                            .rev()
                            .map(|h| h.0)
                            .collect::<Vec<Tau>>();
        cps.push(*t);

        let mut tk1 = cps[1];
        let mut m_tk = mu0.clone();
        let mut params = Vec::with_capacity( (*k as usize)+1 );
        let ct_sigma2 = ChangeType::new("Step", &[sigma02])?;

        // 管理状態
        let ct_mu0 = ChangeType::new("Step", &[mu0])?;
        params.push( (tk1, ct_mu0, ct_sigma2) );

        if *k > 0 {
            for tk in cps[2..].iter() {
                let beta_tk = betak.get_beta_k(&tk1, tk, &m_tk)?;
                let ct_mu_tk = ChangeType::new("Linear", &[beta_tk, m_tk])?;
                params.push( (*tk, ct_mu_tk, ct_sigma2) );
                m_tk = m_tk + beta_tk * ((*tk - tk1) as Mu);
                tk1 = *tk;
            }
        };

        
        // match process_param::norm::Scenario::new(mle.n, &params){
        //     Ok(s) => Ok(s),
        //     Err(e) => Err(e)
        // }
        process_param::norm::Scenario::new(mle.n, &params)
    }



    /// 任意の変化点個数における最適なモデルを取得
    ///
    /// # 引数
    /// * `k` - 変化点個数
    pub fn get_model_with_k(&self, k: &NumChg) -> Result<process_param::norm::Scenario, Box<dyn std::error::Error>> {
        Self::get_model(&self.mle.t, k, &self.incontrol, &self.mle, &self.betak, &self.loglike) 
    }


    /// 任意の変化点個数における最適なモデルを尤度とともに取得
    ///
    /// # 引数
    /// * `k` - 変化点個数
    pub fn get_model_ll_with_k(&self, k: &NumChg) -> Result<(process_param::norm::Scenario, f64), Box<dyn std::error::Error>> {
        let model = Self::get_model(&self.mle.t, k, &self.incontrol, &self.mle, &self.betak, &self.loglike)?;
        let val = self.loglike.get_value(&self.mle.t, &k)?;
        Ok((model, val))
    }


    /// 尤度関数を最大にするモデルを取得
    pub fn get_model_max_ll(&self) -> Result<process_param::norm::Scenario, Box<dyn std::error::Error>> {
        let max_k = <LogLike as CalcDPWithVari<f64, Mu, (MleNorm, Betak, Mu)>>::calc_max_k(&self.mle.t);
        let model_ll = (0..=max_k).map(|k| self.get_model_ll_with_k(&k))
                                  .collect::<Result<Vec<(process_param::norm::Scenario, f64)>, Box<dyn std::error::Error>>>()?;
        let (model_max_ll, _) = model_ll.iter()
                                        .reduce(|acc, ml| {
                                            if acc.1 >= ml.1 {
                                                acc
                                            } else {
                                                ml
                                            }
                                        }).unwrap();
        Ok(model_max_ll.clone())
    }


    // AICを用いた評価関数
    fn calc_aic(loglike: &LogLike, n: &Tau, t: &Tau, k: &NumChg, sigma02: &Sigma2) -> Result<f64, CalcDpError> {
       let ll = loglike.get_value(t, k)?;
       Ok( -(*n as f64) * ll / (2.0 * *sigma02) + 2.0 * (*k as f64) )
    }


    // すべての変化点個数に対してAICを計算
    fn calc_aic_all(loglike: &LogLike, n: &Tau, t: &Tau, sigma02: &Sigma2) -> Result<Vec<(Tau, f64)>, CalcDpError> {
        let max_k = <LogLike as CalcDPWithVari<f64, Mu, (MleNorm, Betak, Mu)>>::calc_max_k(t);
        let aics = (0..=max_k).map(|k| Self::calc_aic(loglike, n, t, &k, sigma02) )
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
                sw.append_row(xlsx::row!["k: Position of change points", "Change point","Beta_k", "Mu at tau_{k-1}"])?;
                for (k, mk) in self.model().parameters().iter().enumerate() {
                    let tk = mk.tau() as f64;
                    match mk.mu() {
                        ChangeType::Step{level} => sw.append_row(xlsx::row![(k as f64).to_cell_value(), tk.to_cell_value(), (), (level as f64).to_cell_value()])?,
                        ChangeType::Linear{grad, init} => match init {
                            Some(init_val) => sw.append_row(xlsx::row![(k as f64).to_cell_value(), tk.to_cell_value(), (grad as f64).to_cell_value(), (init_val as f64).to_cell_value()])?,
                            None => sw.append_row(xlsx::row![(k as f64).to_cell_value(), tk.to_cell_value(), (grad as f64).to_cell_value(), ()])?,
                        },
                        // `ChangeType::StLi`は互換性のために残しています．`Err`にしてもいいかもしれません．
                        ChangeType::StLi{grad, init} => sw.append_row(xlsx::row![(k as f64).to_cell_value(), tk.to_cell_value(), (grad as f64).to_cell_value(), (init as f64).to_cell_value()])?,
                    };
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


    /// 任意の変化点個数でのモデルについて，解析結果をExcelで出力
    ///
    /// # 引数
    /// * `xlsx_path` - 保存先のExcelファイルパス
    /// * `k` - 変化点個数
    pub fn to_excel_k(&self, xlsx_path: &str, k: &NumChg) -> Result<(), Box<dyn std::error::Error>> {
        let mut wb = xlsx::Workbook::create(xlsx_path);
        let model_tk = self.get_model_with_k(k)?;

        let mut sheet_mle = wb.create_sheet("Mu");
        { 
            let model_params_all = model_tk.decomplession()?;
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
                sw.append_row(xlsx::row!["k: Position of change points", "Change point","Beta_k", "Mu at tau_{k-1}"])?;
                for (k, mk) in model_tk.parameters().iter().enumerate() {
                    let tk = mk.tau() as f64;
                    match mk.mu() {
                        ChangeType::Step{level} => sw.append_row(xlsx::row![(k as f64).to_cell_value(), tk.to_cell_value(), (), (level as f64).to_cell_value()])?,
                        ChangeType::Linear{grad, init} => match init {
                            Some(init_val) => sw.append_row(xlsx::row![(k as f64).to_cell_value(), tk.to_cell_value(), (grad as f64).to_cell_value(), (init_val as f64).to_cell_value()])?,
                            None => sw.append_row(xlsx::row![(k as f64).to_cell_value(), tk.to_cell_value(), (grad as f64).to_cell_value(), ()])?,
                        },
                        // `ChangeType::StLi`は互換性のために残しています．`Err`にしてもいいかもしれません．
                        ChangeType::StLi{grad, init} => sw.append_row(xlsx::row![(k as f64).to_cell_value(), tk.to_cell_value(), (grad as f64).to_cell_value(), (init as f64).to_cell_value()])?,
                    };
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


    /// 尤度を最大にするモデルをExcelファイルで出力
    ///
    /// # 引数
    /// * `xlsx_path` - 保存先のExcelファイルパス
    pub fn to_excel_max_ll(&self, xlsx_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let model_max_ll = self.get_model_max_ll()?;
        let k = (model_max_ll.parameters().len() - 1) as NumChg;
        self.to_excel_k(xlsx_path, &k)
    }
}
