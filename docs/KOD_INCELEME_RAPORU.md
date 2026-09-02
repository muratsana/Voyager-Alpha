# Voyager Alpha — Eleştirel Kod İnceleme Raporu

**İncelenen sürüm:** `d5e7804 Initial Voyager Alpha release` (D:\Software Development\Voyager Alpha), 11.545 satır Python, PyQt6 + numpy/scipy/astropy/astroquery/sep/exoplanet-core/pyqtgraph.
**Yöntem:** Tüm `core/` modülleri satır satır okundu; `gui/` modülleri yapı ve akış düzeyinde okundu; iki kaide dokümanındaki ([Q] kaynaklı) kurallarla karşılaştırıldı. Test dosyaları var ama çalıştırılmadı (yalnızca kapsam değerlendirildi).
**Genel yargı:** Mimari iskelet (worker thread'ler, kalibrasyon, kayıt, medyan şablon, SkyBoT, NEA TAP önbelleği, Gaia/VSX ön kontrol) doğru yöne bakıyor ve ortalama bir amatör projeden iyi. Ancak **bilimsel çıktı üretebilecek durumda değil**: her iki modülde de sonucun sayısal değerini (transit derinliği / astrometrik konum ve zaman) sistematik olarak bozan hatalar var. Aşağıdaki liste ciddiyet sırasına göredir. Kod alıntıları dosya:satır ile verilmiştir.

Ciddiyet: **K1** = sonucu yanlış yapar / bilimsel kullanımı engeller; **K2** = hassasiyet/doğruluk kaybı, yanlış pozitif/negatif; **K3** = eksik özellik, kalite/rapor; **K4** = temizlik, performans, test.

---

## A. Her iki modülü etkileyen altyapı sorunları

| # | Ciddiyet | Sorun | Yer | Doğru yaklaşım (kaide dokümanı referansı) |
|---|---|---|---|---|
| A1 | **K1** | **Zaman sistemi yanlış/eksik.** `midpoint_jd = DATE-OBS + EXPTIME/2` yalnızca JD_UTC üretir; BJD_TDB dönüşümü hiç yok (`quality_flags` ile itiraf edilmiş: `time_is_jd_utc_not_bjd_tdb`). `TIMESYS`, `DATE-AVG`, `DATE-END`, `JD`, `MJD-OBS` anahtarları okunmuyor; DATE-OBS'un poz sonu yazıldığı kameralar için kontrol yok. Gözlemevi konumu (`SITELAT/SITELONG/SITEELEV`) hiç okunmuyor. | `core/metadata.py:37-45`, `core/exoplanet.py:323` | Ötegezegen dok. §1.9: astropy ile `t.tdb.jd + light_travel_time(...)`; gözlemevi konumu FITS'ten veya kullanıcı profilinden; `AIRMASS` hesabı. Asteroid dok. §9: ADES `obsTime` orta poz UTC, `rmsTime`. |
| A2 | **K1** | **Fotometri ve tespit, yeniden örneklenmiş (warp edilmiş) karelerde yapılıyor.** `warp_affine(order=1)` bilinear interpolasyon piksel gürültüsünü korele eder, PSF'yi genişletir, tepe değeri düşürür; saturasyon testi için ham kareye ayrıca bakılıyor ama akı/merkez/FWHM hep hizalanmış kareden geliyor. | `core/registration.py:393-416`, `core/exoplanet_worker.py:279-321`, `core/pipeline.py:209-226` | Açıklıkları **ham karede** ölç: hedef/karşılaştırma koordinatlarını her kareye affine/WCS ile taşı (`sensor_coordinates` zaten var), ölçümü ham kalibre karede yap. Warp yalnızca görüntüleme, medyan şablon ve blink için (Siril "shift-only" gerekçesi; asteroid dok. §3.3). |
| A3 | **K2** | **Kalibrasyon denklemi hatalı.** Bias çıkarılmadan dark ölçekleniyor (`dark * t_s/t_d`) → dark bias içeriyorsa bias da ölçeklenir. Master flat, bias/dark düzeltilmiş kabul ediliyor ama doğrulanmıyor; flat için `> 0.05` eşiği kötü piksel maskesi üretmiyor. Master kare **üretimi** (ham bias/dark/flat'ten medyan) hiç yok; kullanıcı dışarıdan master getirmek zorunda. Sıcaklık/binning/gain/filtre eşleşme kontrolü yok. | `core/calibration.py:254-303` | Ötegezegen dok. §3.3 denklemi; §3.5 eşleşme tablosu; master üretimi (medyan/3σ-kırpılmış, N≥10 uyarısı); kötü piksel maskesi (flat <0.5/>1.5, dark >med+5σ) → asteroid `build_defect_mask` ile birleştir. |
| A4 | **K2** | **Plate solve kalitesi ölçülmüyor.** ASTAP `-r 30` kör arama, `-fov` ve `-z` yok, `.ini` dosyasındaki `PLTSOLVD`, artık RMS ve yıldız sayısı okunmuyor; çözüm yalnızca "WCS var mı" ile kabul ediliyor. Asteroid modülünde yalnızca **referans kare** çözülüyor, diğer kareler affine ile taşınıyor (`propagate_wcs_header`). | `core/plate_solver.py:349-407`, `core/pipeline.py:406-420`, `core/exoplanet_worker.py:100-117` | Asteroid dok. §3.2: her kare bağımsız çözüm, kabul RMS ≤ 0.5″ ve ≥ 20 yıldız; `-fov` başlıktan (`XPIXSZ/FOCALLEN`), `-z 2`, `-update`; `.ini` parse; astrometry.net yedek. |
| A5 | **K2** | **Kayıt (registration) toleransları gevşek ve dönüşüm modeli belirsiz.** "Güvenilir" eşiği RMS ≤ 2.0 px ve ≥ 8 yıldız — astrometri için 2 px (çoğu kurulumda 2–4″) çok yüksek. Faz korelasyonu yedeği yalnızca öteleme; alan dönmesi/ölçek olduğunda (meridyen dönüşü, kutup hizasızlığı) sistematik hata. | `core/exoplanet_worker.py:103`, `core/pipeline.py:217`, `core/registration.py` | Kayıt yalnızca görüntüleme/şablon için kullanılacaksa gevşek eşik kabul; **ölçüm** için WCS-tabanlı konum aktarımı (A2). Kalite eşiği: RMS ≤ 0.5 px; ≥ 20 yıldız. |
| A6 | **K3** | **FITS başlığından okunmayan kritik alanlar:** `SITELAT/LONG/ELEV`, `AIRMASS`, `OBJCTRA/DEC`, `RA/DEC`, `TIMESYS`, `SATURATE` (yalnızca ötegezegen), `RDNOISE`/`EGAIN` (okuma gürültüsü hiç kullanılmıyor), `OBSERVER`, `TELESCOP`, `APTDIA`. | `core/metadata.py` | Ötegezegen dok. §1.10 tablosu; asteroid dok. §2.7 "dizi yetenek kartı" için gerekli. |
| A7 | **K3** | **Ağ erişimi senkron ve sessizce yutuluyor.** Gaia/VSX/SkyBoT/NEA istekleri `except Exception: return` ile yutulur; kullanıcı "kontrol yapılamadı"yı görmez, zaman aşımı 12–15 s UI thread'i kilitleyebilir (bazıları worker içinde, bazıları değil). Önbellek yok (SkyBoT `_cache` yalnızca oturum içi). | `core/exoplanet_quality.py:499-568`, `core/known_objects.py:147-159`, `core/known_objects.py:315-334` | Tüm ağ çağrıları worker'da, hata sınıflandırmalı (offline / timeout / HTTP / parse), log'a `WARN` ve sonuç bayrağı; disk önbelleği (SQLite) alan+epoch anahtarlı; kullanıcıya "çevrimdışı mod" göstergesi. |
| A8 | **K3** | **VSX uç noktası eski** (`www.aavso.org/vsx/index.php`) — `vsx.aavso.org`'a taşındı; JSON şeması `VSXObjects.VSXObject` doğru. Gaia sorgusu `TOP 500` ve `launch_job_async` (async gereksiz, sync yeter ve daha hızlı). | `core/exoplanet_quality.py:507-544` | Ötegezegen dok. §4.7. |
| A9 | **K4** | `core/astrometry.py:check_known_asteroid` ölü kod (kullanılmıyor, `known_objects` ile çakışıyor); `detection.build_static_sky_model` kullanılmıyor (memmap sürümü var); `plate_solver._looks_solved` kullanılmıyor. | — | Sil / birleştir. |
| A10 | **K4** | **Test kapsamı** sentetik, çoğunlukla "kod çalışıyor mu" düzeyinde; hiçbir test bilinen doğru sonuca (ör. sentetik geçiş derinliği 1.0%, sentetik asteroid hızı 0.5″/dk) **sayısal tolerans** ile doğrulamıyor; gerçek FITS regresyon verisi yok. | `voyager_alpha/tests/*` | Görev listesi T-Q1..Q4. |

---

## B. Ötegezegen (Exoplanet) modülü

### B1. Detrending geçişi "yiyor" — **K1**

`_robust_polynomial_trend` (`core/exoplanet.py:492-512`) geçiş bilinmeden, **zamana karşı 2. derece polinomu tüm ışık eğrisine** uyduruyor. Alt %30 akı noktalarını atmak (`threshold = percentile(flux, 30)`) ve asimetrik kırpma (−2.5σ / +4σ) geçici bir yama: geçiş süresi verinin %30'undan uzunsa (çok yaygın: 2 saat geçiş, 4 saat gözlem = %50) polinom geçişin içine oturur ve derinliği küçültür/şeklini bozar. Ayrıca detrend **zamana** karşı; fiziksel değişken **hava kütlesi** hiç yok (hesaplanmıyor bile, bkz. A1).

Doğru yaklaşım (ötegezegen dok. §5.8, TFOP SG1 kuralı): (1) OOT bölgesi katalog efemerisinden veya kullanıcı işaretinden tanımlanır; (2) trend yalnızca OOT'ye fit edilir **ya da** geçiş modeliyle **ortak** fit edilir (`χ² = Σ[(O − Σc_j D_j − E)/σ]²`); (3) detrend parametreleri AIRMASS-önce, ΔBIC > 2 kuralıyla en fazla 2 tane; (4) kısmi geçişte detrend kapalı.

### B2. Model fiti fiziksel önsellerden kopuk — **K1/K2**

`fit_limb_darkened_transit` (`core/exoplanet.py:358-452`):
- Limb darkening sabit `(0.3, 0.2)`; filtre ve yıldız Teff/logg'a göre değil (ExoTETHyS/LDTk/Claret; dok. §5.9). Rc bandındaki K cücesi ile B bandındaki F yıldızı aynı katsayıyı alıyor → derinlik yanlılığı.
- Periyot, a/R*, eğim ve eksantriklik yok; "süre" → "hız" dönüşümü `velocity = 2√((1+k)²−b²)/T14` ile geometrik ama `impact` üst sınırı 1.2 (grazing dışı fiziksel değil) ve `T14` ile `T23` ayrımı yok.
- Katalogdan önsel alınmıyor (`pl_ratror, pl_ratdor, pl_orbincl, pl_orbper` NEA'da var, `exoplanet_catalog` zaten indiriyor ama saklamıyor).
- `snr = depth / (median_error/√(n/4))` uydurma; `delta_bic` null modeli 1 parametreli sabit — fakat detrend zaten önceden uygulandığı için BIC karşılaştırması önyargılı.
- Belirsizlik yok: `least_squares` kovaryansı bile kullanılmıyor; MCMC/nested sampling yok → `T_mid ± σ`, `Rp/R* ± σ` rapor edilemiyor (ExoClock/AAVSO gönderimi için zorunlu).
- Hata ölçekleme (`fit_error = max(error, 0.35·median)`) keyfi; ExoClock sırası (ML fit → 3σ kırp → RMS'e ölçekle → β → MCMC) yok.

### B3. Açıklık fotometrisi kaidelere uymuyor — **K2**

`aperture_measurement` (`core/exoplanet.py:108-186`):
- Açıklık yarıçapı FWHM'den bağımsız sabit (`4–10 px` combo). Kaide: 1.5–2 × FWHM taraması, min OOT RMS seçimi (dok. §5.1).
- Halka: `r_in = max(r+3, 1.5r)`, `r_out = max(r_in+4, 2.3r)` — `N_sky` ≥ N_ap garantisi yok (r=6: N_ap≈113, halka 9–14 px ≈ 361 px, tamam; r=10: N_ap≈314, halka 15–23 ≈ 955, tamam; ama küçük r'de sınırda).
- Arka plan medyan+MAD, yıldız maskesi yok, sigma-clip yok (kaide: 3σ kırpılmış medyan, halkadaki tespit edilmiş yıldızlar maskeli).
- Gürültü: `√(F/G + n_pix σ_bg²(1+n_pix/n_sky))` — okuma gürültüsü, dark ve **scintillation** yok (dok. §5.5–5.6 CCD denklemi + Osborn). Gain yalnızca `EGAIN` veya ≤ 10 ise; aksi halde 1 e⁻/ADU varsayımı → hata çubukları yanlış.
- FWHM ikinci moment ile (doğru ama gürültülü); merkezleme kütle merkezi + %8 tepe eşiği ile, 2-B Gauss yok (dok. §5.3).
- Piksel örtüşmesi ikili (`radius <= r`) — `exact`/alt-piksel yok; küçük açıklıklarda %2–5 akı sıçraması.

### B4. Karşılaştırma ensemble'ı ve ön kontrol — **K2** (kısmen iyi)

`exoplanet_quality.py` kaidelerin çoğunu içeriyor (Δmag ±1.5 → ideal −0.44…+0.75, Gaia BP−RP, VSX, leave-one-out, komşu). Eksikler:
- Ağırlık `1/σ²` tek geçiş; Broeg yinelemesi ve mesafe ağırlığı yok (dok. §2.2 Adım 8). Kabul edilebilir ama belgelenmeli.
- Karşılaştırma **adayları hedefin çevresindeki 36 yıldızla** sınırlı ve yalnızca **referans karede** saptanıyor; sürüklenme sonrası kenar dışına çıkanlar kontrol edilmiyor (`edge_margin` yalnızca ilk kare).
- `_apply_leave_one_out_stability` sınırı `max(0.005, 3×medyan)` — mutlak %0.5 alt sınır, iyi kurulumda (RMS 1 mmag) gerçek değişkenleri kaçırır; teorik σ'ya (CCD denklemi) göre `k=2` olmalı (Astrokit).
- Renk farkı eşiği 0.7 → uyarı; kaidede ≤0.3–0.5 tercih, ≤1.0 sınır. Tamam.
- `analysis_allowed` ≥ 2 referans; kaide ≥ 3 (uyarı var).
- Ön kontrol **tüm kareleri iki kez** okuyor (preflight + analiz) → yavaş; ölçüm önbelleği yok.

### B5. Aday karar mantığı ve kapsama — **K2**

`differential_light_curve_from_fluxes` (`core/exoplanet.py:311-322`):
- `coverage_ok = start ≥ 2 and n−end ≥ 2` → yalnızca **2 nokta** taban çizgisi yeterli sayılıyor; kaide ≥ 30 dk / ≥ 0.5 T_dur her iki tarafta (dok. §1.1).
- `depth ≥ max(0.002, 3σ)` + `SNR ≥ 5` + (`ΔBIC ≥ 6` veya box SNR ≥ 7): makul ama ExoClock kabul kriterleri (derinlik S/N ≥ 3, Rp/R* literatürle 3σ içinde, Gauss artık, β) yok.
- Aykırı değer temizliği yok (3σ hareketli medyan; dok. §5.7).
- Kalite metrikleri yok: β (Winn 2008), RMS-vs-bin, χ²_red, ETD DQ; `scatter > 0.02` tek uyarı.
- Box araması `max_width ≤ n/3` → uzun geçişleri (verinin >%33'ü) bulamaz.

### B6. Katalog entegrasyonu — **K3** (iyi temel)

`exoplanet_catalog.py` NEA TAP'tan 4 tablo, SQLite, atomik değişim: **doğru**. Eksikler:
- Sütunlar: `pl_orbpererr1, pl_tranmiderr1, pl_ratror, pl_ratdor, pl_orbincl, pl_orbeccen, sy_gaiamag, pl_refname` yok → efemeris belirsizliği ve fit önselleri üretilemiyor.
- `toi` tablosunda `toidisplay` sütun adı **doğrulanmalı** (NEA TOI sütunları: `toi, toipfx, tid, ...`); yoksa parser sessizce `TOI-` üretir.
- ExoClock `planets_json` (öncelik, güncel efemeris, min teleskop), ExoFOP CSV, exoplanet.eu (Controversial/Retracted) yok (dok. §4).
- Overlay yalnızca **seçili hedef için 5′ koni**; tüm alan overlay'i (WCS footprint) yok. Kullanıcı isteği: "fotoğraftaki bilinen ötegezegenli yıldızlar overlay olarak" → `cone_search(footprint)` + `all_world2pix`.
- `predict_nearest_transit` BJD/JD_UTC karışıklığı (A1 çözülünce düzelir); ±σ_Tmid yok.
- SQLite koni araması dec-bant + tam ayrım taraması; HEALPix/R-tree yok (7k satırda kabul edilebilir, TOI+KOI+K2 ~20k → hâlâ kabul edilebilir).

### B7. Çıktı/rapor — **K3**
- AAVSO Exoplanet Database formatı, ExoClock 3-sütun dosyası, ETD Dmag çıktısı yok; yalnızca CSV.
- Grafik: pyqtgraph tek panel (`Relative flux` vs `Elapsed time`), hata çubuğu var; artık paneli, binli veri, karşılaştırma eğrileri, sistematik paneli (FWHM/sky/airmass/x/y), giriş/çıkış çizgileri, RMS/β anotasyonu yok (dok. §6.1).

---

## C. Asteroid (hareketli nesne) modülü

### C1. Astrometri raporlanamaz durumda — **K1**

- Tespit RA/Dec'i **referans karenin WCS'i** ile, warp edilmiş piksel koordinatından hesaplanıyor (`detection.py:171-177`, `pipeline.py:254`). Kayıt hatası (≤2 px kabul!) doğrudan astrometrik hataya dönüşür; her kare için bağımsız WCS yok (A4).
- ADES taslağında `obsTime = DATE-OBS` (**poz başı**, orta poz değil), `rmsRA/rmsDec/mag/band/astCat/mode/stn` boş; `# telescope`, `# observers` blokları yok (`core/ades.py`). MPC bunu reddeder.
- σ_konum hiç hesaplanmıyor (`FWHM/(2.355·SNR)` + plate-solve RMS; dok. §2.5).
- Konum açısı **piksel uzayında** `atan2(vx, vy)` (`tracklet.py:129`): görüntü yönelimine bağlı, gökyüzü PA (Kuzey→Doğu) değil; SkyBoT `dRA/dDEC` ile karşılaştırılamaz.

### C2. Bilinen nesne eşleme yetersiz — **K2**

`known_objects.py`:
- SkyBoT sorgusu **geosentrik** (`location` verilmiyor → 500); NEO'larda paralaks onlarca ″ hata. MPC gözlemevi kodu veya `lat/lon/alt` verilmeli.
- Yalnızca **referans kare epoğunda** sorgu; tracklet'in tüm zaman aralığı boyunca pozisyon tahmini yok. Eşleme `tolerance_arcsec=20` sabit; SkyBoT `Err` sütunu ve hız/PA uyumu (`|Δhız| ≤ %20, |ΔPA| ≤ 15°`) kullanılmıyor (dok. §7.3).
- `objFilter`, `filter` parametreleri yok; kuyruklu yıldız/gezegen ayrımı `Class` sütunundan yalnızca string.
- JPL sb_ident / Horizons / NEOCP-Scout / MPChecker yedekleri yok; çevrimdışı MPCORB yolu yok.
- "Görünürlük" (`measure_local_peak`) tepe piksel/MAD ile, açıklık fotometrisi değil; `limiting_magnitude` Gaia G ile ölçülüyor (iyi fikir) ama V↔G dönüşümü yok.

### C3. Tracklet bağlama: hız penceresi ve tutarlılık kuralları eksik — **K2**

`tracklet.py`:
- Hız penceresi **piksel/kare** (`min 1.5 px/kare`, `max 35 px/kare`): plate scale ve kadansa göre ″/dk'ya çevrilmiyor; kullanıcı "ana kuşak / NEO / TNO" rejimi seçemez (dok. §2.6, §5 Adım 1).
- Doğrusal fit **ağırlıksız**; χ²_red yok; PA/hız segment tutarlılığı (≤10°/≤%20), kadir tutarlılığı (≤0.5 mag), SNR/FWHM oranı (0.5–2×) yok.
- `max_residuals_per_frame = 24` — yoğun alanda gerçek nesneleri keser; SNR sıralı kesim yerine tüm tespitler + KD-tree kullanılmalı (`scipy.spatial.cKDTree`).
- `_frame_times_minutes` zaman monoton değilse **sessizce kare indeksine düşüyor** (`tracklet.py:325-326`) → yanlış hız; hata vermeli.
- Seed çifti yalnızca 3 kare aralığına kadar (`max_seed_gap_frames=3`); uzun dizide kaçırır.
- `confidence` formülü keyfi; digest2/Find_Orb entegrasyonu yok; NEOCP gönderilebilirlik kriterleri yok.
- Sınıflandırma yalnızca `fit_rms ≤ 0.9 px` → "unknown_candidate": 0.9 px sabit; plate scale'e göre ″ cinsinden olmalı (≤1.0″ MBA, ≤1.5″ NEO).

### C4. Tespit ve şablon — **K2/K3**

`detection.py`:
- Medyan şablon "exact temporal median" — iyi. Ama seeing değişiminde PSF eşleme yok (ZOGY/HOTPANTS seçeneği yok) → her yıldızın çevresinde dipol artıklar; `_has_negative_dipole` bunu bastırmaya çalışıyor (yama).
- `sep.extract(thresh=σ)` sabit `err` skaler; arka plan haritası `sep.Background` ile hesaplanıp `centered` çıkarılmış — tamam. `deblend_cont=0.01` makul.
- Saturasyon testi **artık (residual) görüntüde** (`_looks_saturated`) — anlamsız; ham karede `L_lin` ile yapılmalı.
- Sıcak piksel/kozmik ışın: FWHM alt sınırı + destek pikseli — makul; ama `defect_mask` yalnızca master dark varsa; flat'ten kötü piksel yok.
- Kalıcı artık bastırma (`suppress_persistent_residuals`) hücre 2.5 px, `minimum_hits = max(3, 12%·N)` — yavaş nesneleri (dizi boyunca < 0.7 px hareket) sabit sayar: bu kabul edilebilir bir tasarım kararı ama **kullanıcıya bildirilmeli** (v_min etkisi).
- İz (streak) modu: `ratio ≥ 1.8 & area ≥ 6` sezgisel; Vereš erf-iz fiti yok; iz PA'sı tracklet PA'sı ile karşılaştırılmıyor.

### C5. Sentetik izleme — **K3**

`synthetic_tracking.py`: yalnızca **seçili tracklet çevresinde** (radius 28 px) hız taraması; alan genelinde kör arama yok (Tycho/KBMOD modu). Sabit yıldız maskesi stack **öncesi** uygulanmıyor (kaide §6 Adım 1). Stack SNR = tepe/MAD (radius ≤4 px) — KBMOD Ψ/Φ olabilirliği veya en azından açıklık SNR'ı olmalı. Kare başına doğrulama (≥ %60 karede ≥1.5σ) yok. `synthetic_max_motion_px_hour=12` px/saat sabit.

### C6. Kullanıcı tarafı eksikler
- "Dizi yetenek kartı" (hangi hız/kadir aralığı görülebilir; iz kaybı; okuma gürültüsü rejimi) yok → kullanıcı poz/kadans hatasını sonuç çıkmayınca anlıyor.
- Hız vektörü/hareket ekseni overlay'i yalnızca "seçili tracklet izi"; bilinen nesneler için beklenen yol ve ok yok (kullanıcı isteği).
- Uydu ayrımı (TLE/sat_id) yok; meteor/uçak çizgisi ayrımı yok.
- Blink: hız 1–20 fps var; kayıtlı kareler arası "fark" modu var; ancak tracklet'e kilitli takip (object-centered blink) yok.
- Gözlemevi kodu/konum profili yok (yalnızca ADES'te `UNKNOWN`).

---

## D. Gözlemci (kullanıcı) tarafında eksik olanlar — yazılımın istemesi gerekenler

1. **Gözlemevi profili:** enlem/boylam/yükseklik, MPC kodu (yoksa `XXX` + koordinat), AAVSO gözlemci kodu, saat senkron yöntemi (NTP/GPS) ve son senkron zamanı.
2. **Kamera profili:** doğrusallık sınırı (ölçülmüş ADU; program bir "doğrusallık testi sihirbazı" sunmalı), gain e⁻/ADU, okuma gürültüsü e⁻, dark akımı, full well, piksel boyutu, Bayer düzeni.
3. **Teleskop profili:** açıklık, odak uzaklığı, f/, plate scale doğrulaması (plate solve ile karşılaştırılır), filtre seti (AAVSO/ADES band kodları ile eşleme).
4. **Gözlem planı girdisi:** hedef (katalogdan seç) → beklenen T1/T4, önerilen poz/kadans, taban çizgisi; asteroid için rejim (MBA/NEO/TNO) → hız penceresi, poz sınırı.
5. **Ham kalibrasyon kareleri** (bias/dark/flat setleri) — program master üretmeli; kullanıcıdan master istemek yerine.
6. **Zaman doğrulama kanıtı:** kullanıcının NTP kontrolünü işaretlemesi; FITS zamanlarının düzenliliği otomatik kontrol.

---

## E. Öncelikli düzeltme sırası (özet)

1. A1 zaman sistemi (BJD_TDB, orta poz, airmass, konum) — her şeyin temeli.
2. A2 ham karede fotometri/astrometri; kayıt yalnızca görüntüleme.
3. B1+B2 ortak fit (detrend + Mandel-Agol, LD tablosu, katalog önselleri, MCMC/nested, β).
4. A4+C1 her kare plate solve + σ_konum + doğru ADES.
5. C3 hız penceresi (″/dk), ağırlıklı fit, tutarlılık, KD-tree; C2 SkyBoT topo + hız/PA eşleme + sb_ident.
6. B3 FWHM tabanlı açıklık taraması, CCD denklemi + scintillation, 2-B Gauss merkez.
7. B5/B7 kalite metrikleri ve standart grafik/rapor; ExoClock/AAVSO çıktıları.
8. A3 master üretimi ve kötü piksel maskesi; A7 ağ katmanı.
9. GUI yeniden tasarımı (ayrı doküman).

Her madde `AI_CODER_GOREV_LISTESI.md` içinde kabul kriterli görev olarak açılmıştır.
