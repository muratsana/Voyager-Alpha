# Voyager Alpha

Voyager Alpha, amatör astronomların ardışık FITS gözlemlerini iki ayrı bilimsel
iş akışında incelemesi için geliştirilmiş PyQt6 masaüstü çalışma istasyonudur.

## Modüller

### Asteroid Hunter

1. FITS zaman ve metadata kalite kontrolü
2. İsteğe bağlı master bias, poz süresine göre ölçeklenen master dark ve normalize master flat kalibrasyonu
3. WCS yoksa referans kareyi ASTAP ile otomatik çözme
4. Yıldız eşleştirmeli affine veya alt piksel faz korelasyonlu hizalama
5. Tüm hizalı karelerden disk önbellekli exact temporal-median statik gökyüzü modeli
6. PSF-komşuluk, sıcak piksel ve sensör-sabit kusur filtresinden sonra gerçek zaman aralıklarını kullanan tracklet bağlantısı
7. SkyBoT ile alan/zaman bazlı bilinen cisim kurtarma ve ayrı bilinmeyen cisim taraması
8. Tam kaliteli blink, Synthetic Track doğrulaması, insan incelemesi ve ADES taslak raporu

Bir **tracklet**, aynı hareketli objeye ait olduğu hesaplanan, zaman sıralı çoklu
kare konum ölçümleri zinciridir. Tek karedeki parlak bir kaynak asteroid kanıtı
değildir. `ASTEROID CANDIDATE: YES` yalnızca çoklu kare hareket kanıtı bulunan ve
bilinen obje eşleşmesi olmayan aday üretildiğini söyler; keşif onayı değildir.

### Exoplanet Transit

1. En az beş zaman damgalı FITS karesi
2. Görüntü üzerinden hedef yıldız ve bir veya daha fazla karşılaştırma yıldızı
3. ASTAP referans çözümü, affine WCS propagation ve zayıf kayıtta kare bazlı ASTAP fallback
4. İsteğe bağlı bias/dark/flat kalibrasyonu ve alt piksel kare hizalama
5. Kare bazında yıldız merkezleme, annulus arka planı ve nokta belirsizliği
6. Kararlılığa göre ağırlıklandırılmış çoklu karşılaştırma yıldızı ensemble'ı
7. Robust sabit/doğrusal/ikinci derece detrending ve tek-transit box araması
8. `exoplanet-core` ile quadratic limb-darkened fiziksel transit uyumu
9. Derinlik belirsizliği, SNR, süre, merkez zamanı, `Rp/R*`, impact parameter ve ΔBIC raporu

### Transit Katalogları

Program içi katalog yöneticisi NASA Exoplanet Archive TAP servisinden dört veri
kümesini yerel SQLite cache'e atomik olarak günceller:

- PSCompPars doğrulanmış transit gezegenleri
- TESS TOI / ExoFOP disposition kayıtları
- Kepler KOI confirmed, candidate ve false-positive kayıtları
- K2 planets, candidates, false-positive ve refuted kayıtları

Seçili hedef WCS ile RA/Dec'e çevrilir ve beş yaydakika içindeki katalog kayıtları
`VERIFIED`, `CANDIDATE`, `UNVERIFIED` veya `FALSE POSITIVE` olarak gösterilir. Katalog
ephemeris'i BJD, mevcut fotometri zamanı JD_UTC olduğundan transit pencere karşılaştırması
yaklaşık sonuç olarak etiketlenir. Güncelleme hatasında son başarılı kaynak verisi silinmez.

Transit kararı da aday sınıflandırmasıdır. Airmass, değişken karşılaştırma yıldızı,
meridyen geçişi, bulut ve zaman doğruluğu gibi sistematikler ayrıca incelenmelidir.

## Görüntüleme

- `Auto STF`, `Asinh`, `Manual STF` ve `Linear` ekran germe kipleri
- 1-20 fps önizleme önbellekli, her kareye bağımsız Auto STF uygulayan blink/timelapse
- Seçili tek FITS karesini ASTAP ile plate solve etme
- Renk indeksli overlayler: yeşil bilinen cisim, kırmızı residual aday ve camgöbeği seçili tracklet
- FITS header kaynaklı kamera, sensör, piksel boyutu, binning, gain/offset, sıcaklık, filtre ve optik örnekleme özeti
- Renk kodlu `INFO`, `WARN`, `ERROR` analiz günlüğü
- Tüm temel kontrollerde açıklayıcı mouseover araç ipuçları

## Gereksinimler

- Python 3.10+
- Windows 11
- `numpy`, `scipy`, `exoplanet-core`, `astropy`, `astroquery`, `sep`, `PyQt6`, `pyqtgraph`
- Plate solve için önerilen ASTAP yolu: `C:\Program Files\astap\astap.exe`

## Çalıştırma ve Test

```powershell
.\run_app.ps1
.\run_tests.ps1
.\run_diagnostics.ps1
```

## Windows EXE

Tek dosyalık taşınabilir sürüm:

```powershell
.\build_portable.ps1
```

Çıktı: `dist\Voyager-Alpha.exe`

## Yöntem ve Lisans Sınırı

Bilimsel iş akışı, Citizen Astronomy asteroid/comet kılavuzunda belgelenen Generate,
Align, temporal-median residual, hybrid point/streak Discover, lineer tracklet,
Potential Discovery / Borderline Review, Gaia görünürlük limiti, bilinen cisim
kurtarma ve Synthetic Track kurallarını bağımsız uygular. Varsayılanlar 5 sigma,
3 px FWHM, kare başına 24 residual, 6 px kenar payı, en az 3 bağlı kare, 1.5 px
seed hareketi, 2.8 px yeniden eşleme ve 0.9/1.8 px RMS sınıflarıdır. Kaynak depo
`CC BY-NC-ND` lisanslı olduğundan kodu veya arayüzü kopyalanmamıştır.

Ötegezegen fiziksel model hesabında MIT lisanslı `exoplanet-core` NumPy arayüzü
kullanılır. Tam `exoplanet`/PyMC MCMC yığını masaüstü tarama akışına eklenmemiştir;
bu sürüm hızlı, deterministik aday incelemesine odaklanır.

## Bilimsel Sınır

MPC/ADES gönderimi veya yeni asteroid/ötegezegen iddiası için uygulama çıktısı tek
başına yeterli değildir. Astrometrik artıklar, doğru gözlemevi kodu, zaman sistemi,
fotometrik sistematikler, katalog sorguları ve bağımsız tekrar gözlemleri uzman
kontrolünden geçmelidir. ADES çıktısı bu nedenle taslak olarak işaretlenir.
