# cpd_norm

プロセスが正規分布に従う場合（$ (\bar{X}, s) $ 管理図等が用いられる場合）に対する変化点検出手法。
現状では，次の手法を実装：

- 複数回のステップ変化が生じる場合：`step_multiple`モジュール
- 複数回の線形変化が生じる場合：`linear`モジュール

## 製作者

Author: Shuto Tanabashi / 棚橋秀斗  
Mail: [tanabashi@s.okayama-u.ac.jp](tanabashi@s.okayama-u.ac.jp)  

## 注意事項

本プログラムを実行あるいは流用したことにより生じる一切の事象について製作者は責任を負いません。
自己責任でお願いいたします。

## 実行方法

シェルを開きます。
（WindowsならPowerShellやWindows Terminal。macOSやLinuxならzsh, bash, ターミナル。）

cdコマンド等を使用して`cpd_norm_p`ディレクトリまで移動します。

次のコマンドでこのプロジェクトのドキュメントが読めます。  
注：rustをインストールしていない場合は[rustup](https://www.rust-lang.org/ja/tools/install)をインストールしてください。

```zsh
cargo doc --no-deps --open
```

次のコマンドで実行できます。
実行すると`test/test_scenario.toml`で定義されたシナリオに応じた乱数列が適宜生成され、その乱数列に対して複数回の線形変化を想定した変化点検出が実行されます。
実行結果はターミナルエミュレータ上で表示されるほか、`test`ディレクトリ内に`test_scenario_analyzed_linear.xlsx`として保存されます。

```zsh
cargo run
```
