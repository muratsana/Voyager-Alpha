# Asteroid / Hareketli Nesne Tespiti — Yazılım Kaideleri ve Algoritmalar

**Sürüm:** 1.0 · **Tarih:** 2 Eylül 2026 · **Kapsam:** Aynı alanın art arda çekilmiş karelerinde (CCD/CMOS) hareket eden nokta kaynakların (asteroid, kuyruklu yıldız, NEO, TNO) tespiti, bilinen nesnelerin plate-solve ile overlay edilmesi, aday nesnelerin hareket ekseni/yönü/hızının gösterilmesi ve MPC'ye rapor edilebilir astrometri üretimi. Ötegezegen dökümanının (`otegezegen_gecis_fotometrisi_kaideleri.md`) kardeşidir; kalibrasyon, plate solve ve FITS başlığı kuralları oradan aynen geçerlidir.

**Etiketleme:** **[Q]** = kaynaktan alıntı/aktarım; **[R]** = kaynakların uzlaşısından türetilen spec önerisi. Yalnızca [Q] değerleri "kanıtlanmış", [R] değerleri kullanıcı ayarlanabilir varsayılan olarak kodlanmalıdır.

**Not — meteor değil asteroid:** Meteorlar tek karede çizgi (streak) bırakır ve kataloğu yoktur; bu doküman **kareler arası hareket eden nokta kaynakları** ele alır. Tek karede uzun çizgi bırakan nesneler (uydu, uçak, hızlı NEO) Bölüm 8'de yalnızca **ayırt etme** amacıyla işlenir.

---

## İçindekiler

0. Mimari ve veri akışı
1. Görüntü alma kaideleri (asteroid için farklı olanlar)
2. Enstrüman parametrelerinin etkisi — piksel ölçeği, odak uzaklığı, poz süresi, gürültü, iz kaybı
3. Kare hazırlığı — kalibrasyon, kayıt (registration), plate solve, sabit kaynak eleme
4. Kaynak çıkarma (per-frame detection)
5. Tracklet oluşturma — hareket tutarlılığı ve puanlama (ana algoritma)
6. Shift-and-stack / sentetik izleme (sönük nesne modu)
7. Bilinen nesne kaynakları, overlay ve "bilinmeyen aday" bayrağı
8. Yanlış pozitifler ve uydu/kozmik ışın/sıcak piksel ayrımı
9. Astrometri, fotometri ve zamanlama hassasiyeti
10. Raporlama — ADES, MPC 80-sütun, ALCDEF, NEOCP kuralları
11. Grafik/overlay kaideleri
12. Kabul kontrol listesi
13. Kaynaklar

---

## 0. Mimari ve veri akışı

```
[Kare dizisi: N ≥ 3 (tracklet) / N ≥ 11 (shift-and-stack)]
        │
        ▼
[Kalibrasyon (bias/dark/flat) + kötü piksel maskesi]            ← ötegezegen dokümanı Bölüm 3
        │
        ▼
[Her kare bağımsız plate solve (Gaia DR3); RMS ≤ 0.5″, ≥ 20 yıldız]
        │
        ├──► [Bilinen nesne sorgusu: SkyBoT / JPL sb_ident / MPChecker / NEOCP-Scout] ──► OVERLAY (Bölüm 7, 11)
        │
        ▼
[Sabit kaynak eleme: Gaia maskesi + medyan şablon çıkarma (veya ZOGY/HOTPANTS fark görüntüleme)]
        │
        ▼
[Kare başına kaynak çıkarma: SNR ≥ 5 (rapor) / ≥ 3.5–4 (arama katmanı); sharpness/roundness kesimleri]
        │
        ▼
[Tracklet oluşturma: KD-tree çift → doğrusal fit → 3./4. kare uzatma → artık/hız/PA/kadir tutarlılığı]
        │
        ▼
[Puanlama + bilinen nesne eşleme] ──► [Bilinmeyen aday: Find_Orb Väisälä fit, digest2, sat_id kontrolü]
        │
        ▼
[Astrometri (σ_RA/σ_Dec), fotometri (band, photCat), zaman (mid-exposure UTC)]
        │
        ▼
[Çıktı: ADES PSV/XML, MPC80, ALCDEF, overlay PNG, tracklet tablosu]
Sönük nesne modu: [Sabit kaynak maskesi → (v, PA) ızgarasında shift-and-stack → ≥ 7.5–10σ → kare başına doğrulama]
```

---

## 1. Görüntü alma kaideleri (asteroid için farklı olanlar)

Ötegezegen çekim kuralları (örnekleme ≥ 3 px FWHM, doğrusallık, saat senkronu, FITS başlığı) aynen geçerlidir. Asteroid için **ek** kurallar:

| Konu | Kural | Etiket / kaynak |
|---|---|---|
| Kare sayısı | Tracklet: **≥ 3 kare** ("pratik minimum", yanlış tespitleri önlemek için); yüksek güven 4–5; MPC en az **2 gözlem/nesne/gece** kabul eder ama 3+ tercih | [Q] Petit+ 2004; Astrometrica MOD ≥ 3; MPC Guide |
| Zaman aralığı | Kareler arası ≥ 15 dk (MOPS TTI) veya beklenen hareket ≥ 2 px; toplam yay ≥ 30–45 dk (ana kuşak), **≥ 2 saat** NEOCP'ye gönderilecek NEO için (< 2 s yaylar ZTF'de ~%50 kaybediliyor) | [Q] MOPS; ZTF; [R] |
| Poz süresi (iz kaybı) | `ω · t_poz ≤ FWHM` → iz kaybı < %12; BAA: `t_max(dk) = FWHM(″)/hız(″/dk)` | [Q] Zhai+ 2024; BAA |
| Alt-poz (stacking) | Gökyüzü-sınırlı olmalı: `t_alt ≥ RON²/gökyüzü_oranı(e⁻/px/s)` | [Q] ESO/Hainaut CCD denklemi |
| Sentetik izleme | **≥ 11 kare** minimum, **20+** önerilir; dithering "kesinlikle gerekli" (Catalina ekibi) | [Q] MPC Tycho Tracker sayfası |
| Takip | Yıldız (siderik) takip; nesne hızı biliniyorsa nesne takibi (non-sidereal) yalnızca fotometri/derin stack için | [R] |
| Dithering | Sentetik izleme için gerekli (sıcak piksel/sabit artefakt ayrımı); tracklet modunda isteğe bağlı ama faydalı | [Q]/[R] |
| Zamanlama | Mid-exposure UTC; ana kuşak için ≤ 1 s yeterli; **NEO (1–10°/gün) için ≤ 0.1 s**; hız 3″/s ise 1/6 s hata ≈ 0.5″ konum hatası | [Q] Project Pluto GPS notu |
| Astrometrik katalog | **Gaia DR3** (MPC tercihi; astCat `Gaia3`); USNO-A/B, UCAC1–3, B1950 kabul edilmez; öz hareket görüntü epoğuna taşınmalı | [Q] MPC Guide |
| Plate solve kalitesi | Artık RMS ≤ 0.3–0.5″; her kare bağımsız çözülür (WCS kopyalanmaz) | [Q] EURONEAR 0.3″; Tycho < 0.5″; [R] |
| Alan seçimi | Ekliptik ±15° bandı ana kuşak yoğunluğu; karşı konum (opposition) civarı en parlak; Ay'dan uzak | [R] |

---

## 2. Enstrüman parametrelerinin etkisi

Kullanıcının sorusuna doğrudan yanıt: **evet, piksel boyutu, odak uzaklığı ve gürültü tespiti belirleyici biçimde etkiler.** Yazılım bu parametreleri FITS başlığından (`XPIXSZ`, `FOCALLEN`, `XBINNING`, `GAIN`, `EGAIN`, `RDNOISE`) veya kullanıcı profilinden okuyup aşağıdaki türetilmiş büyüklükleri hesaplamalı ve her gözlem için "bu dizi hangi hız/parlaklık aralığını görebilir" raporu üretmelidir.

### 2.1 Türetilmiş büyüklükler

```
plate_scale [″/px] = 206.265 × piksel_boyutu[µm] × binning / odak_uzaklığı[mm]
FOV [′]            = plate_scale × N_px / 60
FWHM_px            = FWHM_arcsec / plate_scale          → hedef 2–3 px (≥ 3 px astrometri için ideal)
```
- **Örnekleme:** FWHM < 2 px → sıcak piksel/kozmik ışın ile yıldız ayrımı zorlaşır (Astrometrica "Minimum FWHM" 0.7–2 px filtresi ancak gerçek yıldızlar ≥ 2 px ise çalışır) [Q]; FWHM > 5 px → SNR kaybı (n_pix ∝ FWHM²). "2″/px veya daha iyi astrometri için faydalıdır" [Q].
- **Odak uzaklığı** doğrudan plate scale'i ve FOV'u belirler: kısa odak → geniş alan → daha çok asteroid ama düşük astrometrik hassasiyet; uzun odak → tersi. Yazılım FOV'dan beklenen asteroid sayısını (SkyBoT sorgusu) göstererek kullanıcıyı yönlendirir.

### 2.2 Gürültü ve sınır kadiri (CCD denklemi)

```
SNR = S·t / sqrt( S·t + n_pix·(Sky·t + Dark·t + RON²) )        S, Sky: e⁻/s ; n_pix ≈ π·(1.5·FWHM_px)²
σ_mag ≈ 1.0857 / SNR      (SNR 5 → 0.2 mag; SNR 10 → 0.1 mag)
```
- Okuma gürültüsü baskınsa (`Sky·t < RON²`) pozu uzat: `t_min = RON²/Sky` — kısa alt pozlarla stack yapılacaksa bu **zorunlu** kural [Q] Hainaut.
- Sınır kadiri: SNR = 5 (rapor katmanı) veya 3.5–4 (arama katmanı) veren S değerine karşılık gelen kadir; sıfır noktası plate-solve edilen Gaia yıldızlarından (G → V dönüşümüyle) kare başına türetilir.

### 2.3 İz (trailing) kaybı — hız × poz süresi

```
L = ω · t_poz                                (iz uzunluğu; ω ″/s, L ″)
ε_t = SNR_iz / SNR_nokta = w / (w + L)       (w = PSF genişliği ≈ FWHM)   [Q] Shao+ 2014
L ≤ FWHM  →  kayıp < %12                                                     [Q] Zhai+ 2024
```
- Yazılım, seçilen hız penceresinin üst sınırı için `t_poz,max = FWHM/ω_max` hesaplar; kullanıcının pozu bunu aşıyorsa "bu diziyle > ω_kritik hızlı nesneler iz bırakır; iz modunu aç veya alt pozlarla shift-and-stack kullan" uyarısı verir.
- İzli nesneler için nokta açıklığı yerine **Vereš+ 2012 erf-iz modeli** (Gauss ⊗ L uzunluklu çizgi) uydurulur; aspect ≤ 3 için 3×PSF açıklık yeterli, aspect > 10 için iz fiti ~3× daha iyi astrometri ve 2 mag'a kadar daha iyi fotometri verir [Q]. İvme yanlılığı `a·T²/8 ≲ 0.1″` → yakın NEO'larda dikkate alınır [Q].

### 2.4 Stack kazancı

```
SNR_N = sqrt(N) · SNR_1  →  Δm = 1.25 · log10(N)    (11 kare ≈ +1.3 mag, 20 ≈ +1.6, 100 ≈ +2.5)   [Q] Shao 2014 (√N), [R] türetim
```

### 2.5 Astrometrik hassasiyet

```
σ_merkez ≈ FWHM / (2.355 · SNR)   (eksen başına, Gauss PSF)   → FWHM 3″, SNR 10 → 0.13″; SNR 5 → 0.25″   [R]
σ_toplam² = σ_merkez² + σ_plate_solve² + (ω · σ_zaman)²
```
Bu değerler ADES `rmsRA`, `rmsDec` alanlarına yazılır (Bölüm 10).

### 2.6 Hız rejimleri ve birim dönüşümü

```
1 °/gün = 2.5 ″/dk = 150 ″/saat = 0.0417 ″/s
TNO:          1–6 ″/saat            (KBMOD 1–5.7; Petit 1–10)              [Q]
Ana kuşak:    0.2–0.7 °/gün = 0.5–1.75 ″/dk  (MOPS "slow" 0.3–0.7)          [Q]
NEO:          1–5 °/gün (MOPS üst sınır 5; WMOPS 3.25; ZStreak > 4 °/gün iz modu)  [Q]
Uydu (LEO):   > 1 °/saat → neredeyse kesin uydu; GEO ≈ 15 ″/s sabit RA'da    [R]
Minimum hareket: dizi boyunca ≥ 2″ (WMOPS) veya ≥ 2–3 px (Tycho)          [Q]
```

### 2.7 Yazılımın kullanıcıya göstereceği "dizi yetenek kartı" [R]

```
girdi: plate_scale, FWHM, N_kare, t_poz, Δt_toplam, RON, Sky, gain, zero-point
çıktı:
  - Tespit edilebilir en yavaş hız: ω_min = 2 px·plate_scale / Δt_toplam  (″/dk)
  - İzsiz en hızlı hız:             ω_iz  = FWHM_arcsec / t_poz            (″/s → °/gün)
  - Tek kare sınır kadiri (SNR 5) ve N-kare stack sınır kadiri (+1.25 log N)
  - Beklenen astrometrik σ @ SNR 5/10/20
  - Zaman hassasiyeti gereksinimi: σ_t ≤ 0.5″/ω_max
  - Uyarılar: okuma gürültüsü baskın (t < RON²/Sky), undersampled (FWHM < 2 px), kadans > 15 dk vb.
```

---

## 3. Kare hazırlığı

### 3.1 Kalibrasyon
Ötegezegen dokümanı Bölüm 3 aynen. Ek: kötü piksel maskesi ve sıcak piksel haritası burada **daha kritik**; sıcak pikseller dithering yoksa sabit kalır ve sabit kaynak elemeyle gider, dithering varsa "hareket eden" gibi görünebilir → çok kare çakışma testi (Bölüm 8) zorunlu.

### 3.2 Plate solve
Her kare bağımsız (ASTAP/astrometry.net, Gaia DR3 indeksleri); kabul: RMS ≤ 0.5″, ≥ 20 eşleşen yıldız, SIP/3. derece distorsiyon [R]. Öz hareket düzeltmesi görüntü epoğuna (Gaia J2016.0'dan). Kare kalite kapısı: FWHM > 1.5 × medyan veya gökyüzü > 2 × medyan olan kareler stack'ten dışlanır [R].

### 3.3 Kayıt (registration)
Tespit her karenin **kendi piksel uzayında** WCS ile yapılır (interpolasyon gürültüsü yok). Görüntüleme/stack için Siril'in gerekçesiyle **yalnızca kaydırma** (shift-only) dönüşümü tercih; interpolasyon gerekirse Lanczos [Q Siril].

### 3.4 Sabit kaynak eleme (dört yöntem, veri kalitesine göre)

1. **Katalog maskesi:** Gaia DR3 yıldızlarını kare sınır kadirinin ~1 mag altına kadar `r = max(1.5·FWHM, 3 px)` yarıçapla maskele; ~10 mag'dan parlak yıldızlar için kırınım çubuğu/taşma maskesi [R].
2. **Medyan şablon çıkarma:** Siderik kayıtlı N ≥ 5 karenin medyanı (dizi boyunca > 1 FWHM hareket eden nesneler medyanda kaybolur); `kare − medyan`; ≥ (N−1) karede aynı piksel konumunda kalan artıklar sıcak piksel/sabit artık [R]. WMOPS analoğu: 90 dk'da < 2″ hareket → sabit [Q].
3. **Fark görüntüleme:** Seeing dizide > %20 değişiyorsa ZOGY (kapalı form D ve S puan görüntüsü) veya HOTPANTS (`-ng 3 6 0.70 4 1.50 2 3.00 -r 10 -rss 15`) [Q].
4. **Çok kare çakışma:** ≥ 3 karede, farklı zamanlarda, 1 px içinde tekrar eden tespit → sabit → at (sıcak pikselleri de öldürür) [R].

---

## 4. Kaynak çıkarma (per-frame)

```python
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
_, med, std = sigma_clipped_stats(img_masked, sigma=3.0)
f = DAOStarFinder(fwhm=fwhm_px, threshold=5.0*std,          # rapor katmanı 5σ; arama katmanı 3.5–4σ
                  sharplo=0.2, sharphi=1.0, roundlo=-0.5, roundhi=0.5,   # [Q] DAO varsayılan sharp 0.2–1.0; [R] round ±0.5
                  peakmax=0.8*L_sat)
src = f(img_masked - med)
```
Alternatif: `sep.extract(data, thresh=5*bkg.globalrms, minarea=3)`.

**Kesimler [R, kaynak: Astrometrica/MOPS/WMOPS]:**
- SNR ≥ 5 rapor; 3.5–4σ yalnızca ≥ 4 karede doğrulanırsa (MOPS 3σ arşiv, WMOPS 4σ).
- FWHM_kaynak ∈ [0.6, 1.8] × FWHM_yıldız (kare bazlı); Astrometrica: min FWHM 0.7–2 px sıcak piksel reddi [Q].
- sharpness 0.2–1.0 (sıcak piksel/kozmik ışın reddi); |roundness| ≤ 0.5 iz modu kapalıysa.
- tepe < 0.8 × L_sat; kenar payı ≥ kernel yarı genişliği.
- PSF fit χ²_red > 4 → kozmik ışın/artefakt olarak reddet (WMOPS) [Q].
- İz modu (L > 1 FWHM): Vereš erf-iz modeli; iz PA'sı tracklet PA'sı ile 15° içinde uyuşmalı (MOPS "PA aligned") [Q].
- Her tespit için sakla: `x, y, RA, Dec, t_mid(UTC), flux, SNR, FWHM, sharp, round, σ_x, σ_y, tepe, kare_id`.

---

## 5. Tracklet oluşturma — ana algoritma

Kaynak uzlaşısı: Pan-STARRS MOPS (Denneau+ 2013), Kubica+ 2007 KD-tree bağlama, Petit+ 2004, THOR (Moeyens+ 2021), heliolinx, ZTF MODE.

```
girdi: tespit listeleri D_i (kare i, zaman t_i), hız penceresi [v_min, v_max], PA aralığı (isteğe bağlı), toleranslar
0. Koordinatlar: alan merkezine gnomonik projeksiyon (ξ, η) [arcsec]; her karede KD-tree.
1. ÇİFTLER: her d ∈ D_i için, j > i ve Δt = t_j − t_i ≥ Δt_min olan D_j'de
      v_min·Δt ≤ |d_j − d_i| ≤ v_max·Δt   halkasında adayları al.        # Kubica koni araması
   v_min: dizi boyunca ≥ max(2″, 2 px) toplam hareket (WMOPS/Tycho) [Q]
   v_max: rejime göre (TNO 10″/h; MBA ≤ 1°/gün; NEO ≤ 5°/gün MOPS sınırı) [Q]
2. UZATMA: çiftten doğrusal tahminle k > j karelerinde en yakın tespiti al;
      kabul yarıçapı tol = max(1.0″, 2·σ_merkez, 1 px)                        # THOR 1″; Petit 1 FWHM / 3 px; Find_Orb 2″ üst sınır [Q]
3. FIT: ξ(t), η(t) ağırlıklı (1/σ²) doğrusal fit;
      RMS ≤ 1.0″ (MBA/TNO), ≤ 1.5″ (NEO, yay < 30 dk); χ²_red > 4 → ret
      yay > 2 saat ve NEO hızı → kuadratik fit (MOPS "track") + eğrilik anlamlılık testi
4. MİNİMUM TESPİT: 3 (rapor); 2 = yalnızca "çift adayı" (MOPS: 2'li tracklet'lerin ~%10'u gerçek) [Q]; 4–5 yüksek güven (ZMODE 4, WMOPS 5) [Q]
5. TUTARLILIK (yay < 1 saat):
      ardışık segmentler arası hız değişimi ≤ %20, PA değişimi ≤ 10°          [R]
      kadir saçılımı ≤ 0.5 mag (MODE "akı ile bağlama") [Q]; SNR ve FWHM oranları 0.5–2× (MOPS "benzer morfoloji") [Q]
6. TEKİLLEŞTİRME: ≥ 2 tespiti paylaşan tracklet'leri birleştir; uzun/düşük RMS olanı tut (MOPS collapseTracklets) [Q]
7. ALTERNATİF (yoğun alan): (ξ, η, t) uzayında RANSAC/Hough doğru arama; THOR eşdeğeri: hız-projeksiyonlu çerçevede DBSCAN eps ≈ 18″, min 5 tespit [Q]
```

**Puanlama [R]:**
```
puan = w1·f(n_det: 3→0, 4→+1, ≥5→+2)
     + w2·(1 − RMS/σ_merkez_ort)            (normalize artık)
     + w3·(SNR tutarlılığı)  + w4·(1 − kadir_saçılımı/0.5)
     + w5·(1 − Δhız/%20) + w6·(1 − ΔPA/10°)
     + w7·(yay uzunluğu / 60 dk, üst sınır 1)
     + w8·(beklenen karelerde görülme oranı; KBMOD "valid obs" ≥ %80)
ek: digest2 NEO puanı (RA/Dec/zaman hazır olunca; NEOCP için ≥ 65 [Q]); Find_Orb Väisälä RMS (< 0.5″ hedef)
```

**Çıktı (tracklet kaydı):** `id, tespitler[], hız (″/dk), PA (°), RMS (″), n_det, yay (dk), puan, bilinen_eşleşme (isim, ayrım ″, hız/PA uyumu), digest2, durum ∈ {known, unknown_candidate, rejected}`.

---

## 6. Shift-and-stack / sentetik izleme (sönük nesne modu)

Kaynaklar: Shao+ 2014, Zhai+ 2024, KBMOD (Whidden+ 2019), Tycho Tracker.

```
girdi: gökyüzü-sınırlı alt kareler (t_alt ≥ RON²/Sky), N ≥ 11 (20–100 tercih), dithering'li
1. Sabit kaynakları MASKELE (medyan şablon veya Gaia maskesi) — stack'ten ÖNCE; aksi halde yıldız artıkları sırt oluşturur
2. Hız ızgarası: v ∈ [v_min, v_max], adım Δv = 0.5·FWHM / T_toplam (dizi boyunca kayma ≤ yarım PSF)
                 PA ∈ kullanıcı aralığı, adım ΔPA = Δv/v
   → Shao "1 px/kare" ve Zhai 100×100 ızgarasına eşdeğer
3. Her (v, PA) için kareleri kaydır (alt-piksel, shift-only) ve birleştir:
      sigma-clip ortalama, veya KBMOD olabilirlik: L = ΣΨ/√ΣΦ (Ψ = ters-varyans PSF korelasyonu, Φ = PSF etkin alanı)
4. Eşik: ≥ 7.5σ (Zhai: alan başına ~%2 yanlış pozitif) – 10σ (KBMOD lh_level)   [Q]
5. Doğrulama: (i) kare başına damga: nesne karelerin ≥ %60'ında ≥ 1.5σ; (ii) kare akılarında sigmaG [25,75] kırpma (KBMOD);
              (iii) tepe ofseti < 1 px; (iv) kümeleme (KBMOD cluster_eps 20 px) ile tekilleştirme
6. Maliyet: N_kare × N_v × N_PA × N_px → GPU yolu sağlanmalı (Tycho, KBMOD, JPL)
7. Rapor: stack orta zamanı, ADES notes 'K' (stacked), hız/PA, stack SNR, bireysel karelerde görünürlük
```
Tycho kullanıcı deneyimi: "sensitivity" (varsayılan %50), "granularity", min/max hız, PA aralığı, güven puanı ile sıralı liste; örnekte 509 iz → 16 gerçek → kullanıcı doğrulaması **zorunlu** [Q]. MPC uyarısı: "Find_Orb'da yörünge uyması tespitin gerçekliğini DOĞRULAMAZ" [Q].

---

## 7. Bilinen nesne kaynakları, overlay ve bilinmeyen aday bayrağı

### 7.1 Servisler (doğrulanmış)

| Servis | Kullanım | Uç nokta |
|---|---|---|
| **IMCCE SkyBoT** (perturbe efemeris, günlük güncel) | Alan koni araması — ana kuşak için en iyi tek yanıt | `https://vo.imcce.fr/webservices/skybot/skybotconesearch?-ep=<JD/ISO>&-ra=<deg>&-dec=<deg>&-sr=<arcsec\|15m\|0.25d>&-observer=<MPC>&-mime=json&-output=all&-filter=120&-objFilter=111` |
| **JPL sb_ident** | "Bu alanda bu anda hangi küçük cisimler var" — `two-pass=true` ile Horizons kalitesi | `https://ssd-api.jpl.nasa.gov/sb_ident.api?mpc-code=<kod>&obs-time=YYYY-MM-DD_hh:mm:ss&fov-ra-center=..&fov-dec-center=..&fov-ra-hwidth=..&fov-dec-hwidth=..&vmag-lim=21&two-pass=true` |
| **JPL Horizons** | Seçili nesne için mid-exposure zamanlarında referans konum + 3σ belirsizlik | `https://ssd.jpl.nasa.gov/api/horizons.api?COMMAND='DES=…;'&CENTER='<kod>@399'&TLIST=…&QUANTITIES='1,3,9,19,20,23,24,36,37,38'` |
| **MPChecker** | Web formu; toplu kullanım için önerilmez | `https://minorplanetcenter.net/cgi-bin/checkmp.cgi` (radius 5–300′, limit V 20) |
| **NEOCP / Scout** | Onay bekleyen NEO adayları + belirsizlik bulutu | `https://minorplanetcenter.net/iau/NEO/neocp.txt`; `https://ssd-api.jpl.nasa.gov/scout.api?tdes=..&obs-code=..&orbits=true` |
| **PCCP** | Onay bekleyen kuyruklu yıldızlar | `https://minorplanetcenter.net/iau/NEO/pccp.txt` |
| **Sentry** | Risk listesi etiketi | `https://ssd-api.jpl.nasa.gov/sentry.api` |
| **MPC API** | İsim normalizasyonu (`query-identifier`), gözlem çekme (`get-obs`), WAMO durum | `https://data.minorplanetcenter.net/api/…` (GET + JSON gövde; koni araması **yok**) |

SkyBoT çıktı sütunları: `Num, Name, RA, DEC, Class, Mv, Err(″), d(″), dRAcosDEC, dDEC (″/h), Dg, Dh, faz, elongasyon, x,y,z,vx,vy,vz, epoch`. astroquery: `Skybot.cone_search(coo, rad, epoch, location='500', position_error=120)`.

### 7.2 Çevrimdışı yol
- `MPCORB.DAT` / `mpcorb_extended.json.gz` (günlük), `NEAm00…NEAp15.txt` (NEO'lar için ±15 gün epoklu), `CometEls.txt`, Lowell `astorb.dat` (CEU = güncel efemeris belirsizliği ″ → hata çemberi).
- Yayılım: **pyoorb/OpenOrb** (n-cisim, Horizons ile mas düzeyinde uyum) tercih; **Skyfield** iki-cisim (`mpc.load_mpcorb_dataframe`) yalnızca tanımlama için (epoktan ±30 gün içinde < 1–2″, ±100 günde onlarca ″, NEO yakın geçişinde dakikalar) [Q]; **Find_Orb `fo`** kendi ölçümlerinizin mevcut yaya göre artıklarını (O−C) görmek için.
- H-G parlaklık (Bowell 1989): `V = H + 5 log10(rΔ) − 2.5 log10[(1−G)Φ1 + GΦ2]`, `Φ1 = exp(−3.33 tan^0.63(α/2))`, `Φ2 = exp(−1.87 tan^1.22(α/2))`, G bilinmiyorsa 0.15. Kuyruklu yıldız: `m = M1 + 5 log10 Δ + 2.5 K1 log10 r`.

### 7.3 Eşleme ve bilinmeyen aday kuralı [R]

```
for her tracklet T:
    adaylar = bilinen_nesneler(T.orta_zaman, T.orta_konum, r = max(3·σ_efemeris, 10″))
    eşleşme = argmin ayrım; kabul: ayrım ≤ max(3·σ_efemeris, 10″) VE |Δhız| ≤ %20 VE |ΔPA| ≤ 15°
    if eşleşme: T.durum = known; T.isim = MPC query-identifier ile normalize; O−C hesapla (Horizons referans)
    else:
        T.durum = unknown_candidate
        zorunlu ön kontroller: yüksek güven katmanı (≥ 4 tespit) + Find_Orb Väisälä RMS < 0.5″ + digest2 + sat_id (uydu değil) + NEOCP/PCCP'de yok
        → yalnızca bunlar geçerse "NEOCP'ye gönderilebilir" etiketi; aksi halde "ikinci gece gerekli"
```
MPC: hatalı gözlemler kısa yaya "kolayca uyar" → yörünge fiti kanıt değildir [Q].

---

## 8. Yanlış pozitifler ve ayırt etme

| Kaynak | Belirti | Test |
|---|---|---|
| Sıcak piksel | Aynı piksel konumunda tüm karelerde; dithering'de yıldızlara göre "hareket eder" | Kötü piksel maskesi; piksel koordinatında ≥ 3 karede çakışma → at; FWHM < 1 px ve sharpness > 1 |
| Kozmik ışın | Tek kare, keskin, PSF'ye uymaz | Tek karede ortak yok; sharpness > 1; PSF χ²_red > 4 (WMOPS) [Q]; astroscrappy |
| Uydu | Tek karede uzun çizgi veya kareler arası > 1°/saat; GEO ≈ 15″/s | `sat_id` (Project Pluto, TLE: CelesTrak `gp.php?GROUP=active&FORMAT=tle`, ≤ 2 saatte bir); Skyfield ile TLE hızı ↔ ölçülen hız |
| Uydu parlaması (glint) | Tek karede nokta, sonraki karede yok | Tek tespit → çift/tracklet oluşmaz |
| Kırınım çubuğu / hayalet / halo | Parlak yıldızın yakınında, yıldızla sabit ofsetli | Parlak yıldız maskesi (≥ 10 mag), ofset sabitliği testi |
| Yıldız artığı (kötü çıkarma) | Katalog yıldızının 1–2 px yanında, tüm karelerde | Gaia maskesi; sabit konum çakışması |
| Değişken yıldız | Sabit konum, akı değişir | Sabit → at (hareket şartı) |
| Gerçek asteroid ama izli | Uzun, PA tracklet ile uyumlu | İz modu; PA ≤ 15° uyum |
| Meteor | Tek karede çok uzun, parlaklık boyunca değişen çizgi, başka karede yok | Tek kare + uzunluk ≫ ω_max·t_poz → "meteor/uçak" etiketi, tracklet'e alınmaz |

İki bağımsız dedektörün (örn. DAOStarFinder + SEP) kesişimi yanlış pozitifleri 30–200× azaltır, derinlik kaybı yalnızca 0.1–0.3 mag [Q] Petit 2004.

---

## 9. Astrometri, fotometri, zamanlama

- **Konum:** WCS + SIP ile `all_pix2world`; her tespit için `σ_RA = σ_merkez_x · plate_scale`, `σ_Dec` benzer; plate-solve RMS'i karesel toplanır. Gaia öz hareketi uygulanmış referans; astCat `Gaia3`.
- **Zaman:** `t_mid = DATE-OBS + EXPTIME/2` (DATE-OBS obtüratör açılışı; `TIMESYS=UTC` kontrolü; `DATE-AVG` varsa tercih). ADES `obsTime` = `YYYY-MM-DDThh:mm:ss.ssZ`; MPC80 gün kesri 5 ondalık (0.86 s) / hızlı nesnede 6. `rmsTime` (ADES) raporla. Stack'te stack orta zamanı + not `K`.
- **Fotometri:** açıklık 1.5–2 × FWHM; sıfır noktası Gaia G → V (`G − V = −0.01760 − 0.006860(BP−RP) − 0.1732(BP−RP)²`, σ ≈ 0.045) veya APASS/ATLAS; band `Vj`/`Rc`/`G`/`Sr`… (ADES iki harfli formlar; `C`, `L`, `u` yeni gönderimlerde yasak); `photCat`, `photAp` (″), `logSNR`, `rmsMag`. Gece başına nesne başına en az bir kadir **ölçülmeli** (kopyalanmamalı) [Q]. İzli karelerde iz fiti akısı.
- **Işık eğrisi:** 3–5 güneş-renkli karşılaştırma (ATLAS/APASS), ışık-zaman düzeltmesi, birim mesafe ve faz açısına indirgeme (G = 0.15), FALC Fourier periyot araması; ALCDEF dışa aktarım (Bölüm 10).

---

## 10. Raporlama

### 10.1 MPC gözlem kuralları [Q]
- Gece başına nesne başına **≥ 2 gözlem** (3+ tercih, ≥ 20–60 dk yayılmış); tek gözlem/gece tüm partiyi reddettirir.
- Keşif için **ikinci gece** (tercihen yakın gece çifti) → geçici tanımlama.
- Gece başına en az bir ölçülmüş kadir; UTC saat GPS/NTP ile doğrulanmış; iç tutarlılık < 1″.
- Gözlemevi kodu: iki ayrı gecede ≥ 10 numaralı NEA astrometrisi (`XXX` yer tutucu + site koordinatları) sonra form; gezici gözlemci kodu 247/270.
- Kuyruklu yıldız aktivitesi açıkça görülmüyorsa **kuyruklu yıldız iddia etme** (PCCP uyarısı).
- Tycho/sentetik izleme: "marjinal/gürültülü" tespit gönderme.

### 10.2 ADES (PSV) — zorunlu yol
Gönderim yalnızca `https://minorplanetcenter.net/submit_psv` / `submit_xml` (e-posta değil); `submit.xsd` ile doğrulama.

```
# version=2022
# observatory
! mpcCode C10
# submitter
! name Ad Soyad
# observers
! name Ad Soyad
# measurers
! name Ad Soyad
# telescope
! design Newtonian
! aperture 0.30
! detector CMO
! fRatio 4
# software
! astrometry <program> <sürüm>
! objectDetection <program>
! photometry <program>
permID|provID|trkSub|mode|stn|obsTime|ra|dec|rmsRA|rmsDec|rmsTime|astCat|mag|band|photCat|rmsMag|photAp|logSNR|seeing|exp|nStars|notes|remarks
|2024 YR4||CMO|C10|2026-09-02T21:37:12.40Z|148.6712345|16.3838123|0.25|0.22|0.1|Gaia3|19.85|G|Gaia3|0.12|3.0|1.10|2.8|60.0|45||
```
Alanlar: `permID` (numara), `provID` (açık geçici tanım), `trkSub` (bilinmeyen için ≤ 7 alfanümerik kendi tracklet kimliğiniz), `mode` ∈ {CCD, CMO, VID}, `stn`, `obsTime` (Z zorunlu), `ra/dec` derece ≥ 6 ondalık, `rmsRA` (cos δ dahil), `disc` = `*` keşif, `notes` (`K` stacked).

### 10.3 MPC 80-sütun (eski, uyumluluk)
Sütunlar: 1–5 numara, 6–12 paketli geçici tanım, 13 `*`, 14 not1, 15 `C`(CCD)/`B`(CMOS), 16–32 `YYYY MM DD.ddddd`, 33–44 `HH MM SS.dd`, 45–56 `sDD MM SS.d`, 66–70 kadir, 71 band, 72 katalog kodu (Gaia DR3 harfi `CatalogueCodes.html`'den doğrulanmalı), 78–80 gözlemevi kodu. Başlık: `COD, CON, OBS, MEA, TEL, NET, ACK, AC2, COM`.

### 10.4 ALCDEF (ışık eğrisi)
`STARTMETADATA … ENDMETADATA` bloğu (`OBJECTNUMBER, OBJECTNAME, MPCDESIG, CONTACTNAME, OBSERVERS, SESSIONDATE, SESSIONTIME, FILTER, MAGBAND, DIFFERMAGS, STANDARD, DELIMITER, LTCTYPE, LTCAPP, REDUCEDMAGS, PHASE, PABL, PABB, COMPNAME/COMPRA/COMPDEC/COMPMAG`), ardından `DATA=JD|MAG|MAGERR` satırları, `ENDDATA`. Gönderim `https://alcdef.org`.

### 10.5 Gönderim sonrası
WAMO (`https://data.minorplanetcenter.net/api/wamo`, gövde: `["<trkSub> <stn>", …]`) ile durum; kabul edilen verinin artıkları MPC/NEODyS `.rwo` dosyalarında.

---

## 11. Grafik / overlay kaideleri

**Alan görüntüsü overlay katmanları (açılıp kapanabilir):**
1. Bilinen asteroidler (SkyBoT/sb_ident): sembol boyutu `Err`/3σ belirsizlikle orantılı çember (elips: Horizons Q37 SMAA/SMIA/θ); etiket: isim, V, hız (″/dk), PA; sınıf renkleri (MB gri-mavi, NEA turuncu, PHA kırmızı, Comet yeşil, KBO mor, Planet sarı).
2. **Hareket vektörü:** her bilinen ve aday nesne için `[t_ilk, t_son]` aralığında beklenen/ölçülen yol çizgisi; ok başı hareket yönünde; uzunluk = hız × dizi süresi; dizi boyunca konumlar küçük noktalarla (kare başına) işaretlenir; "şu anki kare" konumu dolu.
3. Tracklet'ler: tespit noktaları (kare numarasıyla), fit çizgisi, ± tolerans bandı, RMS; durum rengi (known yeşil, unknown_candidate sarı, rejected kırmızı/soluk).
4. NEOCP/Scout belirsizlik bulutu (orbit örnekleri) yarı saydam noktalar.
5. Sentry/PHA bayrağı (ikon).
6. Uydu geçişleri (TLE tahmini) kesikli mavi.
7. Sabit kaynak maskesi, kötü piksel maskesi (hata ayıklama).
8. Kuzey/Doğu oku, ölçek çubuğu, FOV, epoch, plate-solve RMS/yıldız sayısı.
9. Blink modu: kayıtlı kareler arası geçiş (≥ 3 kare), hız/yön okuyla; kullanıcı onay/ret düğmeleri (Astrometrica MOD, Tycho gibi: **insan doğrulaması zorunlu**).

**Tanı grafikleri:** (a) hız–PA dağılımı (tespit edilen tracklet'ler + bilinen nesnelerin beklenen değerleri), (b) ξ(t), η(t) fit + artıklar (″), (c) tespit SNR ve kadir vs zaman, (d) plate-solve artık haritası (vektör alanı), (e) shift-and-stack (v, PA) olabilirlik haritası, (f) "dizi yetenek kartı" (Bölüm 2.7), (g) O−C (″) tablo/grafik bilinen nesneler için.

---

## 12. Kabul kontrol listesi

| # | Test | Eşik | Sonuç |
|---|---|---|---|
| 1 | Kare sayısı | < 3 | RET (tracklet modu); < 11 → sentetik izleme kapalı |
| 2 | Plate solve | RMS > 0.5″ veya < 20 yıldız | kare UYARI; > 1″ → kare RET |
| 3 | Zaman | `DATE-OBS`/`EXPTIME` yok | RET; NTP doğrulanmamış → UYARI |
| 4 | Toplam yay | < 30 dk (MBA) | UYARI; NEOCP için < 2 saat → "gönderme" |
| 5 | Kare aralığı | beklenen hareket < 2 px | UYARI (v_min yükselt) |
| 6 | İz kaybı | ω_max·t_poz > FWHM | UYARI, iz modu öner |
| 7 | Gürültü rejimi | Sky·t < RON² | UYARI (poz uzat) |
| 8 | Örnekleme | FWHM < 2 px | UYARI |
| 9 | Kalibrasyon | flat/dark yok | UYARI (sabit artıklar artar) |
| 10 | Kare kalitesi | FWHM > 1.5×medyan veya sky > 2×medyan | kare dışla |
| 11 | Tespit SNR | < 5 rapor; < 3.5 her koşulda | tespit dışla |
| 12 | Tracklet artık | > 1″ (MBA) / > 1.5″ (NEO) / > 2″ | RET |
| 13 | Tutarlılık | Δhız > %20, ΔPA > 10°, Δmag > 0.5 | RET |
| 14 | Sabit çakışma | ≥ 3 karede aynı px | sabit → at |
| 15 | Uydu | sat_id eşleşmesi veya > 1°/saat | "uydu" etiketi |
| 16 | Bilinmeyen aday | < 4 tespit veya Väisälä RMS > 0.5″ veya digest2 yok | "gönderilemez" |
| 17 | Fotometri | gece/nesne başına ölçülmüş kadir yok | UYARI (MPC gereksinimi) |
| 18 | ADES | `submit.xsd` doğrulaması başarısız | RET |

---

## 13. Kaynaklar

**Algoritma ve boru hatları:** Denneau+ 2013 Pan-STARRS MOPS https://arxiv.org/abs/1302.7281 · Kubica+ 2007 https://arxiv.org/abs/astro-ph/0703475 · Petit+ 2004 MNRAS 347, 471 · Masci+ ZTF MODE poster https://web.ipac.caltech.edu/staff/fmasci/ztf/Masci_posterPDC2019.pdf · ZStreak https://arxiv.org/abs/1904.09645 · Moeyens+ 2021 THOR https://arxiv.org/abs/2105.01056 · heliolinx https://github.com/heliolinx/heliolinx · HelioLinC RR https://github.com/bengebre/heliolincrr · NEOWISE WMOPS https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/sec4_5.html · Shao+ 2014 https://arxiv.org/abs/1309.3248 · Zhai+ 2024 https://arxiv.org/abs/2401.03255 · KBMOD Whidden+ 2019 https://arxiv.org/abs/1901.02492, parametreler https://epyc.astro.washington.edu/~kbmod/user_manual/search_params.html · Vereš+ 2012 (iz modeli) https://arxiv.org/abs/1209.6106 · ZOGY https://arxiv.org/abs/1601.02655 · HOTPANTS https://github.com/acbecker/hotpants · digest2 https://arxiv.org/abs/1904.09188

**Amatör yazılımlar:** Astrometrica ayarlar http://www.astrometrica.at/Papers/Astrometrica-Settings.pdf · EURONEAR Astrometrica kılavuzu http://www.euronear.org/manuals/Astrometrica-UsersGuide-EURONEAR.pdf · Tycho Tracker https://www.tycho-tracker.com/ · MPC Tycho notu https://docs.minorplanetcenter.net/mpc-ops-docs/observations/tycho-tracker/ · Find_Orb https://www.projectpluto.com/find_orb.htm, environ.def https://github.com/Bill-Gray/find_orb/blob/master/environ.def · GPS zamanlama https://www.projectpluto.com/gps_ast.htm · sat_id https://www.projectpluto.com/sat_id.htm · Siril asteroid https://siril.org/tutorials/asteroid_hunting/ · BAA NEO takip https://britastro.org/asteroids/Astrometry%20of%20NEO's%20-%20Follow-up%20observations.htm · MPO Canopus / MPB 46, 412

**MPC / JPL / IMCCE:** MPC Guide to Minor Body Astrometry https://minorplanetcenter.net/iau/info/Astrometry.html · ADES https://minorplanetcenter.net/iau/info/ADES.html, alan değerleri https://www.minorplanetcenter.net/iau/info/ADESFieldValues.html, ADES-Master https://github.com/IAU-ADES/ADES-Master · 80-sütun https://www.minorplanetcenter.net/iau/info/OpticalObs.html · MPC API https://minorplanetcenter.net/mpcops/documentation/ · MPChecker https://minorplanetcenter.net/cgi-bin/checkmp.cgi · NEOCP https://minorplanetcenter.net/iau/NEO/toconfirm_tabular.html · MPCORB https://www.minorplanetcenter.net/iau/MPCORB.html, format https://www.minorplanetcenter.net/iau/info/MPOrbitFormat.html · SkyBoT https://ssp.imcce.fr/webservices/skybot/api/conesearch/ · astroquery IMCCE https://astroquery.readthedocs.io/en/latest/imcce/imcce.html · JPL sb_ident https://ssd-api.jpl.nasa.gov/doc/sb_ident.html · Horizons https://ssd-api.jpl.nasa.gov/doc/horizons.html · Scout https://ssd-api.jpl.nasa.gov/doc/scout.html · Sentry https://ssd-api.jpl.nasa.gov/doc/sentry.html · SBDB https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html · Lowell astorb https://asteroid.lowell.edu/astorb/ · OpenOrb https://github.com/oorb/oorb · Skyfield https://rhodesmill.org/skyfield/kepler-orbits.html · CelesTrak https://celestrak.org/NORAD/documentation/gp-data-formats.php · ALCDEF https://alcdef.org/docs/ALCDEF_Standard.pdf · Gaia DR3 fotometrik ilişkiler https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu5pho/cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html

**CCD fiziği:** Hainaut CCD SNR https://www.eso.org/~ohainaut/ccd/sn.html · Dhillon CCD denklemi https://vikdhillon.staff.shef.ac.uk/teaching/phy217/detectors/phy217_det_ccdeqn.html · photutils DAOStarFinder https://photutils.readthedocs.io/en/stable/api/photutils.detection.DAOStarFinder.html

*Boşluklar:* Astrometrica'nın MOD iç toleransları ve Tycho'nun kesin hız/PA/puan varsayılanları erişilebilir sayfalarda yayımlanmamıştır → kullanıcı ayarlı, yukarıdaki [R] varsayılanlarıyla. MPC `data` sayfası ve NEOCP efemeris CGI form alanları bu oturumda doğrudan doğrulanamadı (403/hata) → kodlamadan önce kontrol edilmeli.
