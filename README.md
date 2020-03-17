# lapack-test

## openblas

blasを最適化したもの

### openblasのインストール

```sh
sudo apt install libopenblas-base
sudo apt install libopenblas-dev
```

### blasのライブラリの選択

```sh
sudo update-alternatives --config libblas.so.3
```

### プログラムのコンパイル

```sh
cd test
gcc test-openblas.c -L/usr/include/openblas-base -I/usr/lib/openblas-base -lopenblas -lpthread -lrt
```