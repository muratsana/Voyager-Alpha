# Voyager Alpha — AI Coder Görev Listesi (Issue Formatı)

Bu dosya doğrudan bir AI kodlama ajanına (Claude Code, Codex, Cursor…) verilmek üzere hazırlanmıştır. Her görev: **Amaç · Dosyalar · Yapılacaklar · Kabul kriterleri · Referans** (kaide dokümanı bölümü). Öncelik: P0 (bilimsel geçerlilik), P1 (doğruluk/hassasiyet), P2 (özellik/rapor), P3 (kalite/temizlik). Görevler bağımlılık sırasına göre listelenmiştir; P0'lar bitmeden P2'lere geçilmemelidir.

**Ajan için genel talimatlar:**
- Kaynak: `D:\Software Development\Voyager Alpha\voyager_alpha`. Python 3.12, PyQt6 6.7, numpy/scipy/astropy/astroquery/sep/pyqtgraph. Yeni bağımlılıklar: `photutils`, `batman-package` (veya mevcut `exoplanet-core`), `ldtk` **veya** `exotethys`, `emcee`, `dynesty`, `astroscrappy`, `pyoorb` (isteğe bağlı), `skyfield`.
- Her görevde önce **birim testi** yaz (sentetik veri, sayısal tolerans), sonra kodu değiştir; `run_tests.ps1` yeşil kalmalı.
- Fonksiyon imzalarını koruyarak yeni davranışı **bayrakla** ekle (`use_raw_frame_photometry=True` gibi), GUI'yi kırma.
- Her sayısal varsayılan için kod yorumunda kaide referansı bırak: `# ref: otegezegen_kaideleri §5.1 (AIJ 1.7×FWHM)`.
- Kullanıcıya görünen metinler Türkçe, `self.tr()` ile.

---

## P0 — Bilimsel geçerlilik

### T-01 Zaman sistemi: orta poz UTC, BJD_TDB, hava kütlesi, gözlemevi konumu
- **Dosyalar:** `core/metadata.py`, `core/models.py`, yeni `core/timeframe.py`, `core/exoplanet.py`, `core/ades.py`
- **Yapılacaklar:** `DATE-OBS` + `EXPTIME` (ve `DATE-AVG`, `DATE-END`, `JD`, `MJD-OBS`, `TIMESYS` kontrolü) → `t_mid_utc`; `SITELAT/SITELONG/SITEELEV` (yoksa profil) → `EarthLocation`; hedef koordinatı (WCS merkezi veya seçili yıldız) → `bjd_tdb = t.tdb.jd + light_travel_time(kind='barycentric')`; `AIRMASS` yoksa `AltAz` ile hesapla; `FrameRecord`'a `t_mid_utc, jd_utc, bjd_tdb, airmass, alt_deg, sun_alt_deg` ekle. Zaman monoton değilse **hata** (sessiz indeks fallback'i kaldır: `tracklet._frame_times_minutes`).
- **Kabul:** Test: bilinen bir (RA, Dec, UTC) için `bjd_tdb − jd_utc` astropy referansı ile < 1 ms; `time_is_jd_utc_not_bjd_tdb` bayrağı kalkar; ışık eğrisi `time_system="BJD_TDB"`.
- **Ref:** ötegezegen §1.9, asteroid §9.

### T-02 Ham karede ölçüm; warp yalnızca görüntüleme
- **Dosyalar:** `core/exoplanet_worker.py`, `core/exoplanet_quality.py`, `core/pipeline.py`, `core/detection.py`
- **Yapılacaklar:** Fotometri için her karede `sensor_coordinates()` ile koordinatı taşı, `aperture_measurement(raw_calibrated, …)` çağır; asteroid tespitlerinde artık görüntüdeki (x,y)'yi sensöre taşıyıp **ham kalibre karede** yeniden merkezle (2-B Gauss) ve akı/tepe/FWHM'i oradan al; RA/Dec'i **o karenin kendi WCS'i** ile hesapla (T-04).
- **Kabul:** Sentetik dizi (bilinen akı, alt-piksel kaymalar) → warp'lı ölçüme göre RMS ≥ %20 düşer; tepe ADU ham ile aynı.
- **Ref:** ötegezegen §5.3, asteroid §3.3.

### T-03 Detrend + geçiş modeli ortak fit, LD tablosu, katalog önselleri, MCMC, β
- **Dosyalar:** `core/exoplanet.py` (yeniden yaz: `transit_model.py`, `detrend.py`, `quality_metrics.py`), `core/exoplanet_catalog.py`
- **Yapılacaklar:** (1) OOT bölgesi efemeris/kullanıcıdan; `_robust_polynomial_trend`'i kaldır. (2) Model: Mandel-Agol (batman veya exoplanet-core) parametreleri `T_mid, Rp/R*, a/R*, i(b), F0, c_airmass, c_t(, c_t²)`; P, e, ω katalogdan sabit; poz süresi supersampling. (3) LD: filtre + Teff/logg/[Fe/H] → ExoTETHyS/LDTk (çevrimdışı tablo yedeği: Claret kuadratik). (4) Detrend seçimi: AIRMASS önce; parametre ekleme ΔBIC > 2; en fazla 2. (5) Sıra: ML fit → 3σ kırp → hataları artık RMS'e ölçekle → β (10–30 dk bin) → emcee (≥ 5000 adım, R̂ < 1.05) veya dynesty. (6) Çıktı: tüm parametreler ± σ, χ², χ²_red, BIC, RMS ham/binli, β, derinlik S/N, ExoClock kabul bayrakları.
- **Kabul:** Sentetik geçiş (derinlik 1.00 %, hava kütlesi trendi %1) → derinlik 1.00 ± 0.05 %, T_mid hata < 1 dk; eski kodda derinlik önyargısı > %15 olan senaryoda yeni kod < %3.
- **Ref:** ötegezegen §5.8–5.10.

### T-04 Her kare plate solve + kalite + σ_konum + doğru ADES
- **Dosyalar:** `core/plate_solver.py`, `core/wcs_cache.py`, `core/pipeline.py`, `core/ades.py`, yeni `core/astrometry_quality.py`
- **Yapılacaklar:** ASTAP çağrısı: `-fov` (XPIXSZ/FOCALLEN/NAXIS'ten), `-z 2`, `-r` (WCS ipucu varsa 2°, yoksa 30°), `-update`, `-log`; `.ini` parse (`PLTSOLVD`, `CRVAL`, hata metni); çözüm sonrası Gaia DR3 ile artık RMS ve yıldız sayısı ölçümü (cone + `all_world2pix`); kabul RMS ≤ 0.5″, ≥ 20 yıldız (ayarlanabilir); astrometry.net (yerel `solve-field` veya `astrometry` PyPI) yedek. Her kare bağımsız; `propagate_wcs_header` yalnızca yedek ve `method="propagated"` etiketli. `σ_konum² = (FWHM/(2.355·SNR))² + RMS_solve²`. ADES: `obsTime` orta poz UTC `…Z`, `ra/dec` 7 ondalık, `rmsRA` (cosδ dahil) `rmsDec`, `rmsTime`, `astCat=Gaia3`, `mode=CCD/CMO`, `stn`, `mag/band/photCat/photAp/logSNR/seeing/exp/nStars`, `trkSub`, `# telescope/# observers/# measurers/# software` blokları; `submit.xsd` doğrulaması (ADES-Master).
- **Kabul:** Sentetik alan (Gaia yıldızlarıyla) → çözüm RMS < 0.3″; ADES çıktısı ADES-Master `valid` aracından geçer; MPC80 eşdeğeri üretilir.
- **Ref:** asteroid §3.2, §9, §10.

### T-05 Tracklet bağlama: fiziksel hız penceresi, KD-tree, ağırlıklı fit, tutarlılık, PA
- **Dosyalar:** `core/tracklet.py`, `core/discovery_method.py`, `core/models.py`
- **Yapılacaklar:** Hız penceresi ″/dk (rejim: MBA 0.3–1.75, NEO ≤ 12.5, TNO 0.017–0.1; Özel); plate scale ile px'e çevir; `max_residuals_per_frame` kaldır (SNR ≥ 5 tümü), `cKDTree` ile çift arama; ağırlıklı (1/σ²) doğrusal fit ξ(t), η(t) (gnomonik, ″), RMS ≤ 1.0″ (MBA/TNO) / 1.5″ (NEO), χ²_red ≤ 4; segment tutarlılığı Δhız ≤ %20, ΔPA ≤ 10°; kadir saçılımı ≤ 0.5 mag; SNR/FWHM oranı 0.5–2×; min tespit 3 (2 = "çift adayı" ayrı liste); ≥ 4 yüksek güven; sky PA (K→D) WCS ile; puanlama (asteroid §5); seed gap sınırını kaldır (KD-tree ile tüm çiftler, Δt ≥ Δt_min).
- **Kabul:** Sentetik dizi: 5 enjekte nesne (0.4–5 ″/dk) 30 kare, 2000 sahte tespit → ≥ 4/5 bulunur, yanlış pozitif ≤ 1; PA gerçek ile ≤ 5°.
- **Ref:** asteroid §5.

### T-06 Bilinen nesne eşleme: topo SkyBoT + sb_ident + hız/PA + belirsizlik + hareket vektörü
- **Dosyalar:** `core/known_objects.py`, yeni `core/ephemeris_sources.py`, `gui/viewer.py`
- **Yapılacaklar:** SkyBoT `-observer=<MPC kodu>` veya `lat/lon/alt`; `-filter`, `-objFilter`; her tracklet'in **orta zamanı** için sorgu (önbellekli); eşleme ≤ max(3·Err, 10″) VE Δhız ≤ %20 VE ΔPA ≤ 15°; JPL `sb_ident` (two-pass) yedek; seçili nesne için Horizons `TLIST` referans + 3σ elips; NEOCP/Scout listesi; disk önbelleği (SQLite: alan, epoch, kaynak). Overlay: beklenen yol `[t_ilk, t_son]`, ok başı, kare noktaları, belirsizlik elipsi, sınıf renkleri.
- **Kabul:** Bilinen bir asteroid içeren gerçek dizi (kullanıcı verisi) → eşleşme ayrımı < 5″, hız/PA uyumu raporlanır; NEO senaryosunda geosentrik vs topo fark testi.
- **Ref:** asteroid §7, §11.

## P1 — Doğruluk / hassasiyet

### T-07 Açıklık fotometrisi: FWHM taraması, 3σ medyan gökyüzü, 2-B Gauss merkez, CCD denklemi + scintillation, exact örtüşme
- **Dosyalar:** `core/exoplanet.py` → `core/photometry.py` (photutils tabanlı)
- **Yapılacaklar:** `CircularAperture` + `ApertureStats(sigma_clip=3, maxiters=10)` medyan, halkadaki tespitli yıldızlar maskeli; r_ap ∈ [0.8, 3.0]×FWHM 12 adım, `r_in ≥ max(r+2, 1.9·FWHM)`, `N_sky ≥ 3·N_ap`; min OOT RMS seçimi + "daha küçük aynı derinlik" kuralı; değişken açıklık seçeneği; 2-B Gauss merkez (yedek CoM), kare RET: kayma > 0.5·FWHM veya FWHM Δ > %50; hata: AIJ B1 (gain, RON, dark, n_pix/n_b) + Osborn scintillation (D, t, X, h, C_Y=1.5); `method='exact'`.
- **Kabul:** Sentetik yıldız (bilinen akı, Poisson+RON) → akı yanlılığı < %0.5, σ tahmini gerçek saçılıma ±%15.
- **Ref:** ötegezegen §5.1–5.6.

### T-08 Karşılaştırma seçimi: kaide algoritması, Broeg ağırlığı, teorik σ tabanlı eleme
- **Dosyalar:** `core/exoplanet_quality.py`
- **Yapılacaklar:** Aday havuzu tüm alan (kenar payı `r_out+30`), sıkı 0.5–1.5× / gevşek |Δm| ≤ 2; Gaia BP−RP ≤ 0.5 tercih / 1.0 sınır, `ruwe > 1.4` uyarı; VSX yeni uç nokta (`vsx.aavso.org`); ampirik eleme `rms > 2·σ_teorik` (mutlak %0.5 sınırı kaldır); leave-one-out %3 kuralı; Broeg yinelemeli ağırlık (1e-5) + mesafe ağırlığı (isteğe bağlı); puan tablosu ve ret nedenleri GUI'ye.
- **Kabul:** Enjekte değişken karşılaştırma (0.5 % genlik) elenir; ensemble RMS, tek en iyi karşılaştırmadan düşük.
- **Ref:** ötegezegen §2.2.

### T-09 Kalibrasyon: master üretimi, doğru denklem, eşleşme kontrolü, kötü piksel maskesi, kozmik ışın
- **Dosyalar:** `core/calibration.py`, `core/detection.py::build_defect_mask`
- **Yapılacaklar:** `build_master_bias/dark/flat(files, method=median|sigclip)`; denklem §3.3 (bias-çıkarılmış ölçekli dark, flat normalize, 0/NaN→1); eşleşme (binning/gain/readmode RET, ΔT > 2 °C uyarı, filtre RET); BPM = dark > med+5σ | flat < 0.5 | flat > 1.5 → interpolasyon + bayrak; `astroscrappy` isteğe bağlı; kalibrasyon durumu FITS/rapor.
- **Kabul:** Sentetik bias/dark/flat ile düz alan → kalibre kare düzlüğü < %0.3; N < 5 uyarısı.
- **Ref:** ötegezegen §3.

### T-10 Tespit: ham karede saturasyon, PSF eşleme seçeneği, iz fiti, iki dedektör kesişimi
- **Dosyalar:** `core/detection.py`
- **Yapılacaklar:** `_looks_saturated` ham kare + `L_lin`; seeing değişimi > %20 ise HOTPANTS-benzeri kernel eşleme (basit: Gauss konvolüsyonla FWHM eşitleme) seçeneği; iz modu Vereš erf-modeli ile PA ve uzunluk; DAOStarFinder + SEP kesişimi seçeneği; sıcak piksel/kozmik ışın kuralları DAO sharpness/roundness ile.
- **Kabul:** Enjekte izli nesne (L = 3 FWHM) PA ± 10°; kozmik ışın enjeksiyonu %95 elenir.
- **Ref:** asteroid §4, §8.

### T-11 Sentetik izleme (shift-and-stack) tam mod
- **Dosyalar:** `core/synthetic_tracking.py`
- **Yapılacaklar:** Alan geneli (v, PA) ızgarası (Δv = 0.5·FWHM/T, ΔPA = Δv/v), sabit kaynak maskesi stack **öncesi**, sigma-clip ortalama (KBMOD Ψ/Φ isteğe bağlı), eşik ≥ 7.5σ, kare başına doğrulama (≥ %60 karede ≥ 1.5σ), kümeleme; ilerleme + iptal; numpy vektörize, GPU (cupy) isteğe bağlı; ADES notes `K`.
- **Kabul:** SNR_tek = 2 nesne, 20 kare → stack'te ≥ 7.5σ ve doğru hız ±1 adım.
- **Ref:** asteroid §6.

## P2 — Özellik / rapor

### T-12 Katalog genişletme: NEA ek sütunlar, ExoClock, ExoFOP, exoplanet.eu, tüm alan overlay
- **Dosyalar:** `core/exoplanet_catalog.py`, `gui/exoplanet_workspace.py`
- **Yapılacaklar:** NEA `pscomppars/ps` sütunları (`pl_orbpererr1, pl_tranmiderr1, pl_ratror, pl_ratdor, pl_orbincl, pl_orbeccen, sy_gaiamag, pl_refname`); ExoClock `planets_json` (öncelik, güncel efemeris, min teleskop, O−C); ExoFOP TOI/CTOI CSV; exoplanet.eu TAP (Controversial/Retracted); birleşik tablo + çakışma kuralı (§4.6); `toidisplay` sütun adını doğrula; footprint tabanlı koni + katmanlı overlay; σ_Tmid = √(σ_T0² + (n·σ_P)²); "efemeris eski" uyarısı.
- **Kabul:** WASP-12 alanı → overlay'de ev sahibi yıldız işaretli, T1/T4 ± σ hesaplı; Retracted nesne kırmızı/gizli.
- **Ref:** ötegezegen §4.

### T-13 Kalite metrikleri ve kabul kontrol listeleri (her iki modül)
- **Yapılacaklar:** ötegezegen §7 ve asteroid §12 tablolarını `core/acceptance.py` kural motoru olarak uygula (RET/UYARI/bilgi); GUI'de adım şeridi renklenmesi; rapora ekle. ETD DQ, ExoClock kabul testleri, SG1 NEB (isteğe bağlı).

### T-14 Grafik seti (ötegezegen §6.1–6.2, asteroid §11)
- **Yapılacaklar:** pyqtgraph çok panelli: ışık eğrisi (ham/binli/model/T1-T4/meridyen), artıklar, karşılaştırmalar, sistematikler, RMS-vs-bin, corner; asteroid: hız–PA, ξ/η fit, SNR–zaman, plate-solve artık alanı, stack olabilirlik; PNG/PDF dışa aktarım (matplotlib) SG1 başlık/altbaşlık/legend/anotasyon kutusu ile.

### T-15 Dışa aktarım formatları
- **Yapılacaklar:** AAVSO Exoplanet Database dosyası (§6.4), ExoClock 3-sütun + info (§6.5), ETD Dmag (§6.6), AIJ uyumlu ölçüm tablosu (§6.3); asteroid ADES/MPC80 (T-04), ALCDEF (§10.4), HTML raporu genişletme.

### T-16 Dizi yetenek kartı ve çekim uygunluk kartı
- **Yapılacaklar:** asteroid §2.7 (ω_min, ω_iz, sınır kadiri, σ_astro, σ_t gereksinimi, uyarılar) ve ötegezegen §1.11 S/N tahmini + §1.2 poz/kadans/verim kontrolleri → GUI adım 1'de kart; FITS'ten okunamayan değerler profil diyaloğuna yönlendirir.

### T-17 Uydu/meteor ayrımı ve NEOCP kontrolleri
- **Yapılacaklar:** CelesTrak GP TLE (`GROUP=active`, ≤ 2 saatte bir önbellek) + skyfield ile alan/zaman geçiş tahmini; tracklet hızı > 1°/saat → uydu etiketi; tek kare uzun çizgi → meteor/uçak etiketi; NEOCP/PCCP `txt` listesi ve Scout ephemeris; digest2 (isteğe bağlı, `digest2` ikili) ve Find_Orb `fo` Väisälä RMS entegrasyonu; "NEOCP'ye gönderilebilir" kuralı (asteroid §7.3).

### T-18 GUI yeniden tasarımı
- **Yapılacaklar:** `GUI_YENIDEN_TASARIM.md` §3–§7: tokens + QSS, dock iskeleti, adım şeridi, adım widget'ları, ImageViewer katman yöneticisi + hareket vektörü + imleç okuma + nesne-merkezli blink, sonuç kartı, UnitSpinBox, profil yöneticisi (gözlemevi/kamera/teleskop/filtre/araçlar), boş durum/hata kartları, `QSettings`, kısayollar, `tr()` + `.ts`, pytest-qt testleri, offscreen ekran görüntüsü regresyonu.

## P3 — Kalite / temizlik

### T-19 Ağ katmanı
- `core/net.py`: tek `fetch(url, cache_key, ttl)` (SQLite önbellek, zaman aşımı, yeniden deneme, hata sınıfı); tüm servisler bunun üzerinden; çevrimdışı mod; kullanıcıya görünür "son güncelleme" ve hata nedeni.

### T-20 Ölü kod ve tutarlılık
- `core/astrometry.check_known_asteroid`, `detection.build_static_sky_model`, `plate_solver._looks_solved` kaldır; `exoplanet_quality._usable_gain` ile `exoplanet_worker._usable_gain` tekilleştir; `_robust_sigma` üç kopya → `core/stats.py`.

### T-21 Test altyapısı
- Sentetik veri üreteci (`tests/synth.py`): yıldız alanı + Gaia benzeri katalog + enjekte geçiş / hareketli nesne / kozmik ışın / uydu izi; sayısal tolerans testleri (T-02…T-11 kabul kriterleri); gerçek küçük FITS regresyon seti (kullanıcıdan 3 dizi); CI (`run_tests.ps1` + GitHub Actions Windows).

### T-22 Performans
- Ön kontrol ve analizde kareleri bir kez oku (ölçüm önbelleği); memmap stack zaten var; fotometri döngüsünü vektörize; ilerleme/iptal her worker'da; 100 kare × 16 MP dizide analiz < 3 dk hedefi (CPU).

---

## Bağımlılık grafiği (özet)
T-01 → T-03, T-04, T-06, T-12 · T-02 → T-07, T-10 · T-04 → T-05, T-06 · T-07 → T-08, T-03 · T-09 → T-10 · T-13/T-14/T-15 → T-18 (GUI son aşamada tüm çıktıları bağlar).
